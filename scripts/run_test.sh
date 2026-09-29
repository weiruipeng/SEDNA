#!/usr/bin/env bash
# run_test.sh — end-to-end smoke test for SEDNA using the shipped example data.
#
# Inputs (resolved relative to the SEDNA package root, i.e. the directory that
# contains bin/, scripts/, envs/, examples/):
#   examples/acc_list.txt
#   examples/species.txt
#   examples/ERR10878159_clean_downsample_10pct.fq.gz
#
# Outputs are written under a single temporary work directory whose path is
# printed on exit so the user can inspect them.
#
# Usage:
#   bash scripts/run_test.sh [TAXDB_DIR]
#
# If TAXDB_DIR is provided, the script reuses an existing taxonomy database
# (must contain names.dmp, nodes.dmp, and the uncompressed
# nucl_gb.accession2taxid).  If omitted, the script downloads taxdmp itself into
# $TEST_ROOT/taxdmp (requires network access).

set -euo pipefail

log()  { printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Required command not found on PATH: $1"; }
assert_file() { [[ -s "$1" ]] || die "Expected non-empty file is missing or empty: $1"; }

# ----------------------------------------------------------------------------
# 0. Locate the SEDNA package root
# ----------------------------------------------------------------------------
SELF="$(readlink -f "$0" 2>/dev/null \
        || python3 -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' "$0")"
SEDNA_ROOT="$(cd "$(dirname "$SELF")/.." && pwd)"
BIN_DIR="${SEDNA_ROOT}/bin"
EXAMPLES="${SEDNA_ROOT}/examples"

EXAMPLE_ACC="${EXAMPLES}/acc_list.txt"
EXAMPLE_SPECIES="${EXAMPLES}/species.txt"
EXAMPLE_FQ="${EXAMPLES}/ERR10878159_clean_downsample_10pct.fq.gz"

assert_file "$EXAMPLE_ACC"
assert_file "$EXAMPLE_SPECIES"
assert_file "$EXAMPLE_FQ"

# Optional positional arg: pre-existing taxdmp directory
TAXDB_DIR_USER="${1:-}"

# ----------------------------------------------------------------------------
# 1. Test workspace (unique per run)
# ----------------------------------------------------------------------------
STAMP="$(date '+%Y%m%d_%H%M%S')"
TEST_ROOT="${SEDNA_ROOT}/test_${STAMP}"
mkdir -p "${TEST_ROOT}"

REF_DIR="${TEST_ROOT}/ref_cave"
RUN_DIR="${TEST_ROOT}/run_sample"
DAMAGE_DIR="${TEST_ROOT}/bamdam"
LOG="${TEST_ROOT}/test.log"

log "SEDNA package root: ${SEDNA_ROOT}" | tee -a "$LOG"
log "Test work directory:  ${TEST_ROOT}" | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 2. Prepare PATH and permissions (non-destructive: only +x when needed)
# ----------------------------------------------------------------------------
for f in "$BIN_DIR"/sedna "$BIN_DIR"/sedna-*; do
    [[ -x "$f" ]] || chmod +x "$f"
done
export PATH="${BIN_DIR}:$PATH"
log "Added ${BIN_DIR} to PATH" | tee -a "$LOG"

# Verify the dispatcher and key deps
need sedna
need bowtie2
need samtools
need ngsLCA
need bamdam
need python3

sedna --version | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 3. Taxonomy database (reuse existing or download)
# ----------------------------------------------------------------------------
if [[ -n "$TAXDB_DIR_USER" ]]; then
    TAXDB_DIR="$(cd "$TAXDB_DIR_USER" && pwd)"
    log "Reusing existing taxdmp: ${TAXDB_DIR}" | tee -a "$LOG"
    assert_file "${TAXDB_DIR}/names.dmp"
    assert_file "${TAXDB_DIR}/nodes.dmp"
    assert_file "${TAXDB_DIR}/nucl_gb.accession2taxid"
else
    TAXDB_DIR="${TEST_ROOT}/taxdmp"
    mkdir -p "${TAXDB_DIR}" && cd "${TAXDB_DIR}"

    log "Downloading NCBI taxonomy dump into ${TAXDB_DIR}..." | tee -a "$LOG"
    need curl

    curl -fsSL -O https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz names.dmp nodes.dmp
    rm -f taxdump.tar.gz

    curl -fsSL -O \
      https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    gunzip -c nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
    rm -f nucl_gb.accession2taxid.gz

    log "taxdmp download done" | tee -a "$LOG"
fi

# ----------------------------------------------------------------------------
# 4. Stage 1 — build reference database
# ----------------------------------------------------------------------------
log "Stage 1/3: sedna build-db" | tee -a "$LOG"

N_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
log "Using threads=${N_THREADS} (override via OMP_NUM_THREADS or SLURM_CPUS_PER_TASK)" \
    | tee -a "$LOG"

sedna build-db \
    -a "$EXAMPLE_ACC" \
    -o "$REF_DIR" \
    --species-file "$EXAMPLE_SPECIES" \
    --threads "$N_THREADS" \
    2>&1 | tee -a "$LOG" || {
        die "Stage 1 (build-db) failed.  Log tail: ${LOG}"
    }

assert_file "${REF_DIR}/work_build_db/combined.dedup.fasta"
# Bowtie2 index: either .bt2 or .bt2l prefix accepted
if ! compgen -G "${REF_DIR}/ncbi_reference.*.bt2" >/dev/null \
   && ! compgen -G "${REF_DIR}/ncbi_reference.*.bt2l" >/dev/null; then
    die "Bowtie2 index not produced under ${REF_DIR}"
fi
log "Stage 1 OK: merged FASTA + Bowtie2 index present" | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 5. Stage 2 — main alignment + classification workflow
# ----------------------------------------------------------------------------
log "Stage 2/3: sedna run" | tee -a "$LOG"

sedna run \
    -i "$EXAMPLE_FQ" \
    -o "$RUN_DIR" \
    --ref-dir "$REF_DIR" \
    --taxdb-dir "$TAXDB_DIR" \
    -t "$N_THREADS" \
    2>&1 | tee -a "$LOG" || {
        die "Stage 2 (run) failed.  Log tail: ${LOG}"
    }

assert_file "${RUN_DIR}/out.name_sorted.bam"
# LCA: accept any filename the pipeline writes
LCA_FOUND=""
for f in "${RUN_DIR}"/lca_results*; do
    if [[ -s "$f" ]]; then LCA_FOUND="$f"; break; fi
done
[[ -n "$LCA_FOUND" ]] || die "No lca_results* file produced in ${RUN_DIR}"
log "Stage 2 OK: BAM + LCA present (LCA=${LCA_FOUND##*/})" | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 6. Stage 3 — damage analysis with bamdam (species rank, no BAM extract)
# ----------------------------------------------------------------------------
log "Stage 3/3: sedna damage" | tee -a "$LOG"

sedna damage \
    --sedna-dir "$RUN_DIR" \
    --ranks "species" \
    --stranded ds \
    --outdir "$DAMAGE_DIR" \
    2>&1 | tee -a "$LOG" || {
        die "Stage 3 (damage) failed.  Log tail: ${LOG}"
    }

COMPUTE_DIR="${DAMAGE_DIR}/species/compute"
assert_file "${COMPUTE_DIR}/bamdam.species.upto_species.tsv"
assert_file "${COMPUTE_DIR}/bamdam.species.upto_species.subs"
assert_file "${DAMAGE_DIR}/species/read_counts.species.tsv"

N_PLOTS=$(find "${DAMAGE_DIR}/species/damage_plots" -maxdepth 1 -name '*.png' 2>/dev/null | wc -l | tr -d ' ')
log "Stage 3 OK: bamdam TSV + subs + read_counts.tsv present; ${N_PLOTS} damage plot(s)" \
    | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 7. Final summary
# ----------------------------------------------------------------------------
{
    echo
    echo "======================================================="
    echo "SEDNA smoke test PASSED  ($(date '+%F %T'))"
    echo "======================================================="
    echo "Test root     : ${TEST_ROOT}"
    echo "Threads used  : ${N_THREADS}"
    echo "Full log      : ${LOG}"
    echo
    echo "Key artifacts:"
    echo "  ref FASTA   : ${REF_DIR}/work_build_db/combined.dedup.fasta"
    echo "  BAM         : ${RUN_DIR}/out.name_sorted.bam"
    echo "  LCA         : ${LCA_FOUND}"
    echo "  bamdam TSV  : ${COMPUTE_DIR}/bamdam.species.upto_species.tsv"
    echo "  read counts : ${DAMAGE_DIR}/species/read_counts.species.tsv"
    echo "  damage plot : ${DAMAGE_DIR}/species/damage_plots/*.png (${N_PLOTS} files)"
    echo
    echo "To remove test artifacts run:"
    echo "  rm -rf '${TEST_ROOT}'"
    echo "======================================================="
} | tee -a "$LOG"
