#!/usr/bin/env bash
set -euo pipefail

usage() {
cat <<'USAGE'
Usage:
  bash scripts/run_pydamage_by_taxon.sh \
    --sedna-dir output136_downlod \
    --ref-fasta ncbi_reference_download/work_build_db/combined.dedup.fasta \
    --rank species \
    --top 20 \
    --min-reads 20 \
    --threads 16 \
    --outdir sample136.pydamage_by_species

Required:
  --sedna-dir DIR
      Output directory from run_sedna.sh. It must contain:
        out.name_sorted.bam
        abundance.<rank>.tsv

  --ref-fasta FILE
      Reference FASTA used during database construction, for example:
        ncbi_reference_download/work_build_db/combined.dedup.fasta

Optional:
  --rank STR
      Taxonomic rank used to select taxa for plotting:
        species/genus/family/order/class/phylum/superkingdom
      default: species

  --top INT
      Select the top N taxa by read count
      default: 20

  --min-reads INT
      Process only taxa with reads >= this threshold
      default: 20

  --threads INT
      Number of threads
      default: 8

  --case-sensitive
      Make seqkit grep case-sensitive.
      Default: case-insensitive.

  --no-md
      Do not run samtools calmd to add MD tags.
      By default, a coordinate-sorted BAM with MD tags is generated for pydamage.

  --outdir DIR
      Output directory
      default: pydamage_by_taxon

Example:
  bash scripts/run_pydamage_by_taxon.sh \
    --sedna-dir output136_downlod \
    --ref-fasta ncbi_reference_download/work_build_db/combined.dedup.fasta \
    --rank species \
    --top 20 \
    --min-reads 20 \
    --threads 16 \
    --outdir sample136.pydamage_by_species

USAGE
}

log() {
    echo "[$(date '+%F %T')] $*" >&2
}

die() {
    echo "[$(date '+%F %T')] ERROR: $*" >&2
    exit 1
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

SEDNA_DIR=""
REF_FASTA=""
RANK="species"
TOP_N=20
MIN_READS=20
THREADS=8
CASE_SENSITIVE=0
ADD_MD=1
OUTDIR="pydamage_by_taxon"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sedna-dir)
            SEDNA_DIR="$2"
            shift 2
            ;;
        --ref-fasta)
            REF_FASTA="$2"
            shift 2
            ;;
        --rank)
            RANK="$2"
            shift 2
            ;;
        --top)
            TOP_N="$2"
            shift 2
            ;;
        --min-reads)
            MIN_READS="$2"
            shift 2
            ;;
        --threads|-t)
            THREADS="$2"
            shift 2
            ;;
        --case-sensitive)
            CASE_SENSITIVE=1
            shift
            ;;
        --no-md)
            ADD_MD=0
            shift
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown argument: $1"
            ;;
    esac
done

[[ -n "$SEDNA_DIR" ]] || die "Missing required argument: --sedna-dir"
[[ -n "$REF_FASTA" ]] || die "Missing required argument: --ref-fasta"

case "$RANK" in
    species|genus|family|order|class|phylum|superkingdom) ;;
    *) die "Unsupported --rank value: $RANK" ;;
esac

need_cmd python3
need_cmd samtools
need_cmd seqkit
need_cmd pydamage

SEDNA_DIR="$(readlink -f "$SEDNA_DIR")"
REF_FASTA="$(readlink -f "$REF_FASTA")"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

BAM="${SEDNA_DIR}/out.name_sorted.bam"
ABUNDANCE="${SEDNA_DIR}/abundance.${RANK}.tsv"

[[ -s "$BAM" ]] || die "BAM not found: $BAM"
[[ -s "$ABUNDANCE" ]] || die "Abundance table not found: $ABUNDANCE"
[[ -s "$REF_FASTA" ]] || die "Reference FASTA not found: $REF_FASTA"

WORK_DIR="${OUTDIR}/work"
SELECTED_DIR="${OUTDIR}/selected_fasta_${RANK}"
BED_DIR="${OUTDIR}/beds_${RANK}"
BAM_DIR="${OUTDIR}/bam_${RANK}"
PYDAMAGE_DIR="${OUTDIR}/pydamage_${RANK}"
PLOT_DIR="${OUTDIR}/plots_${RANK}"
LOG_DIR="${OUTDIR}/logs"

mkdir -p \
    "$WORK_DIR" \
    "$SELECTED_DIR" \
    "$BED_DIR" \
    "$BAM_DIR" \
    "$PYDAMAGE_DIR" \
    "$PLOT_DIR" \
    "$LOG_DIR"

SELECTED_TAXA="${OUTDIR}/selected_${RANK}_taxa.tsv"

log "SEDNA_DIR: $SEDNA_DIR"
log "BAM: $BAM"
log "REF_FASTA: $REF_FASTA"
log "ABUNDANCE: $ABUNDANCE"
log "RANK: $RANK"
log "TOP_N: $TOP_N"
log "MIN_READS: $MIN_READS"
log "THREADS: $THREADS"
log "OUTDIR: $OUTDIR"

###############################################################################
# Step 1. Prepare the reference FASTA index
###############################################################################

log "[1/6] Checking the reference FASTA index"

if [[ ! -s "${REF_FASTA}.fai" ]]; then
    log "Running samtools faidx"
    samtools faidx "$REF_FASTA"
else
    log "Reference FASTA index already exists; skipping: ${REF_FASTA}.fai"
fi

###############################################################################
# Step 2. Prepare the coordinate-sorted BAM / MD BAM
###############################################################################

log "[2/6] Preparing the BAM for pydamage"

COORD_BAM="${WORK_DIR}/out.coord_sorted.bam"
MD_BAM="${WORK_DIR}/out.coord_sorted.MD.bam"

if [[ ! -s "$COORD_BAM" ]]; then
    log "Converting the name-sorted BAM to coordinate sort order"
    samtools sort \
        -@ "$THREADS" \
        -o "$COORD_BAM" \
        "$BAM" \
        > "${LOG_DIR}/samtools.sort.stdout.log" \
        2> "${LOG_DIR}/samtools.sort.stderr.log"
else
    log "Coordinate-sorted BAM already exists; skipping: $COORD_BAM"
fi

if [[ "$ADD_MD" -eq 1 ]]; then
    if [[ ! -s "$MD_BAM" ]]; then
        log "Running samtools calmd to add MD tags"
        samtools calmd \
            -@ "$THREADS" \
            -b \
            "$COORD_BAM" \
            "$REF_FASTA" \
            > "$MD_BAM" \
            2> "${LOG_DIR}/samtools.calmd.stderr.log"
    else
        log "MD BAM already exists; skipping: $MD_BAM"
    fi

    if [[ ! -s "${MD_BAM}.bai" ]]; then
        samtools index -@ "$THREADS" "$MD_BAM"
    fi

    INPUT_BAM_FOR_EXTRACTION="$MD_BAM"
else
    if [[ ! -s "${COORD_BAM}.bai" ]]; then
        samtools index -@ "$THREADS" "$COORD_BAM"
    fi

    INPUT_BAM_FOR_EXTRACTION="$COORD_BAM"
fi

log "BAM used for extraction: $INPUT_BAM_FOR_EXTRACTION"

###############################################################################
# Step 3. Select taxa from the abundance table
###############################################################################

log "[3/6] Selecting taxa to plot from abundance.${RANK}.tsv"

python3 - "$ABUNDANCE" "$SELECTED_TAXA" "$TOP_N" "$MIN_READS" <<'PY'
import sys
import csv

abundance = sys.argv[1]
out_tsv = sys.argv[2]
top_n = int(sys.argv[3])
min_reads = int(sys.argv[4])

rows = []

with open(abundance, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    required = {"taxid", "name", "rank", "reads", "relative_abundance"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"abundance table missing columns: {missing}")

    for row in reader:
        taxid = row["taxid"].strip()
        name = row["name"].strip()
        rank = row["rank"].strip()

        if taxid in {"unclassified", "NA", "", "0"}:
            continue

        if name in {"unclassified", "NA", ""}:
            continue

        try:
            reads = int(float(row["reads"]))
        except Exception:
            continue

        try:
            rel = float(row["relative_abundance"])
        except Exception:
            rel = 0.0

        if reads >= min_reads:
            rows.append((taxid, name, rank, reads, rel))

rows.sort(key=lambda x: x[3], reverse=True)
rows = rows[:top_n]

with open(out_tsv, "w") as out:
    out.write("taxid\tname\trank\treads\trelative_abundance\n")
    for taxid, name, rank, reads, rel in rows:
        out.write(f"{taxid}\t{name}\t{rank}\t{reads}\t{rel}\n")

print(f"Selected taxa: {len(rows)}", file=sys.stderr)
PY

N_SELECTED=$(tail -n +2 "$SELECTED_TAXA" | wc -l | awk '{print $1}')

if [[ "$N_SELECTED" -eq 0 ]]; then
    die "No taxa were selected. Lower --min-reads or try --rank genus/family."
fi

log "Selected taxa: $N_SELECTED"
log "selected taxa:"
cat "$SELECTED_TAXA" >&2

###############################################################################
# Function: generate a BED file from FASTA
###############################################################################

fasta_to_bed() {
    local fasta="$1"
    local bed="$2"

    python3 - "$fasta" "$bed" <<'PY'
import sys

fasta = sys.argv[1]
bed = sys.argv[2]

records = []

name = None
seq_len = 0

def flush():
    global name, seq_len
    if name is not None:
        records.append((name, seq_len))

with open(fasta, "r", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            # The BAM reference name is usually the ID before the first space in the FASTA header
            name = header.split()[0]
            seq_len = 0
        else:
            seq_len += len(line.strip())

flush()

with open(bed, "w") as out:
    for ref, length in records:
        if length > 0:
            out.write(f"{ref}\t0\t{length}\n")
PY
}

###############################################################################
# Function: run pydamage analyze
###############################################################################

run_pydamage_analyze() {
    local bam="$1"
    local outdir="$2"
    local log_prefix="$3"

    mkdir -p "$outdir"

    set +e
    pydamage analyze \
        --outdir "$outdir" \
        -p "$THREADS" \
        "$bam" \
        > "${log_prefix}.pydamage.analyze.stdout.log" \
        2> "${log_prefix}.pydamage.analyze.stderr.log"
    local status=$?
    set -e

    if [[ "$status" -ne 0 ]]; then
        log "pydamage analyze --outdir failed; trying the -o form"

        set +e
        pydamage analyze \
            -o "$outdir" \
            -p "$THREADS" \
            "$bam" \
            > "${log_prefix}.pydamage.analyze.stdout.log" \
            2> "${log_prefix}.pydamage.analyze.stderr.log"
        status=$?
        set -e
    fi

    return "$status"
}

###############################################################################
# Function: run pydamage plot
###############################################################################

run_pydamage_plot() {
    local csv="$1"
    local outdir="$2"
    local log_prefix="$3"

    mkdir -p "$outdir"

    set +e
    pydamage plot \
        "$csv" \
        --outdir "$outdir" \
        > "${log_prefix}.pydamage.plot.stdout.log" \
        2> "${log_prefix}.pydamage.plot.stderr.log"
    local status=$?
    set -e

    if [[ "$status" -ne 0 ]]; then
        log "pydamage plot --outdir failed; trying the -o form"

        set +e
        pydamage plot \
            "$csv" \
            -o "$outdir" \
            > "${log_prefix}.pydamage.plot.stdout.log" \
            2> "${log_prefix}.pydamage.plot.stderr.log"
        status=$?
        set -e
    fi

    return "$status"
}

###############################################################################
# Step 4-6. For each taxon:
#   extract reference contigs with seqkit
#   extract BAM reads with samtools
#   pydamage analyze + plot
###############################################################################

log "[4/6] Extracting reference contigs by taxon name with seqkit"
log "[5/6] Extracting corresponding reads from BAM"
log "[6/6] Running pydamage analyze + plot"

tail -n +2 "$SELECTED_TAXA" | while IFS=$'\t' read -r TAXID TAXNAME TAXRANK READS REL; do

    SAFE_NAME="$(echo "$TAXNAME" | sed 's/[^A-Za-z0-9._-]/_/g' | sed 's/_\+/_/g' | cut -c1-120)"

    TAXON_FASTA="${SELECTED_DIR}/${TAXID}_${SAFE_NAME}.fasta"
    TAXON_BED="${BED_DIR}/${TAXID}_${SAFE_NAME}.bed"
    TAXON_BAM="${BAM_DIR}/${TAXID}_${SAFE_NAME}.bam"

    TAXON_PYDAMAGE_DIR="${PYDAMAGE_DIR}/${TAXID}_${SAFE_NAME}"
    TAXON_PLOT_DIR="${PLOT_DIR}/${TAXID}_${SAFE_NAME}"

    LOG_PREFIX="${LOG_DIR}/${TAXID}_${SAFE_NAME}"

    log "================================================================"
    log "Processing taxon:"
    log "  taxid: $TAXID"
    log "  name : $TAXNAME"
    log "  rank : $TAXRANK"
    log "  reads: $READS"
    log "  rel  : $REL"

    ###########################################################################
    # 4.1 Build the seqkit grep pattern
    ###########################################################################

    PATTERN_FILE="${WORK_DIR}/${TAXID}_${SAFE_NAME}.pattern.txt"

    python3 - "$TAXNAME" "$PATTERN_FILE" <<'PY'
import sys
import re

name = sys.argv[1]
out = sys.argv[2]

# Escape species/genus names for regex safety.
# This prevents parentheses, dots, and similar characters from being misinterpreted.
pattern = re.escape(name)

with open(out, "w") as f:
    f.write(pattern + "\n")
PY

    ###########################################################################
    # 4.2 Use seqkit to extract reference sequences whose headers contain the taxon name
    ###########################################################################

    SEQKIT_CMD=(
        seqkit grep
        -n
        -r
        -f "$PATTERN_FILE"
        "$REF_FASTA"
    )

    if [[ "$CASE_SENSITIVE" -eq 0 ]]; then
        SEQKIT_CMD+=(-i)
    fi

    log "seqkit command: ${SEQKIT_CMD[*]}"

    "${SEQKIT_CMD[@]}" \
        > "$TAXON_FASTA" \
        2> "${LOG_PREFIX}.seqkit.stderr.log" || {
            log "WARNING: seqkit grep failed: taxid=${TAXID}, name=${TAXNAME}"
            continue
        }

    if [[ ! -s "$TAXON_FASTA" ]]; then
        log "WARNING: no contigs were extracted for this taxon from the reference FASTA:"
        log "  taxid=${TAXID}"
        log "  name=${TAXNAME}"
        log "  Possible reason: the reference FASTA headers do not contain this species or taxon name."
        continue
    fi

    N_CONTIGS=$(grep -c '^>' "$TAXON_FASTA" || true)

    log "Extracted contig count: $N_CONTIGS"
    log "taxon FASTA: $TAXON_FASTA"

    ###########################################################################
    # 4.3 Generate BED
    ###########################################################################

    fasta_to_bed "$TAXON_FASTA" "$TAXON_BED"

    if [[ ! -s "$TAXON_BED" ]]; then
        log "WARNING: BED is empty; skipping: $TAXON_BED"
        continue
    fi

    log "BED: $TAXON_BED"

    ###########################################################################
    # 5. Extract reads from BAM that map to these contigs
    ###########################################################################

    if [[ ! -s "$TAXON_BAM" ]]; then
        log "Extracting reads from BAM"

        # -F 2308 filters:
        #   4    unmapped
        #   256  secondary alignment
        #   2048 supplementary alignment
        # total = 2308
        samtools view \
            -@ "$THREADS" \
            -b \
            -F 2308 \
            -L "$TAXON_BED" \
            "$INPUT_BAM_FOR_EXTRACTION" \
            > "$TAXON_BAM" \
            2> "${LOG_PREFIX}.samtools.view.stderr.log"

        samtools index -@ "$THREADS" "$TAXON_BAM"
    else
        log "Taxon BAM already exists; skipping extraction: $TAXON_BAM"
    fi

    N_ALIGNMENTS=$(samtools view -c "$TAXON_BAM" || echo 0)

    log "Extracted alignment count: $N_ALIGNMENTS"

    if [[ "$N_ALIGNMENTS" -eq 0 ]]; then
        log "WARNING: no alignments were extracted for this taxon; skipping pydamage."
        continue
    fi

    ###########################################################################
    # 6. pydamage analyze
    ###########################################################################

    log "Running pydamage analyze"

    if run_pydamage_analyze "$TAXON_BAM" "$TAXON_PYDAMAGE_DIR" "$LOG_PREFIX"; then
        log "pydamage analyze completed: $TAXON_PYDAMAGE_DIR"
    else
        log "WARNING: pydamage analyze failed:"
        log "  taxid=${TAXID}"
        log "  name=${TAXNAME}"
        log "  Check logs:"
        log "    ${LOG_PREFIX}.pydamage.analyze.stderr.log"
        continue
    fi

    ###########################################################################
    # 6.1 Locate pydamage_results.csv
    ###########################################################################

    PYDAMAGE_CSV=""

    if [[ -s "${TAXON_PYDAMAGE_DIR}/pydamage_results.csv" ]]; then
        PYDAMAGE_CSV="${TAXON_PYDAMAGE_DIR}/pydamage_results.csv"
    elif [[ -s "${TAXON_PYDAMAGE_DIR}/pydamage_results/pydamage_results.csv" ]]; then
        PYDAMAGE_CSV="${TAXON_PYDAMAGE_DIR}/pydamage_results/pydamage_results.csv"
    else
        PYDAMAGE_CSV="$(find "$TAXON_PYDAMAGE_DIR" -type f -name "pydamage_results.csv" | head -n 1 || true)"
    fi

    if [[ -z "$PYDAMAGE_CSV" || ! -s "$PYDAMAGE_CSV" ]]; then
        log "WARNING: pydamage_results.csv was not found; skipping plot:"
        log "  taxid=${TAXID}"
        log "  name=${TAXNAME}"
        continue
    fi

    log "pydamage CSV: $PYDAMAGE_CSV"

    ###########################################################################
    # 6.2 pydamage plot
    ###########################################################################

    log "Running pydamage plot"

    if run_pydamage_plot "$PYDAMAGE_CSV" "$TAXON_PLOT_DIR" "$LOG_PREFIX"; then
        log "pydamage plot completed: $TAXON_PLOT_DIR"
    else
        log "WARNING: pydamage plot failed:"
        log "  taxid=${TAXID}"
        log "  name=${TAXNAME}"
        log "  Check logs:"
        log "    ${LOG_PREFIX}.pydamage.plot.stderr.log"
        continue
    fi

done

log "ALL DONE"
log "Main outputs:"
log "  selected taxa:"
log "    $SELECTED_TAXA"
log "  Taxon FASTA extracted by seqkit:"
log "    $SELECTED_DIR"
log "  BED:"
log "    $BED_DIR"
log "  taxon BAM:"
log "    $BAM_DIR"
log "  pydamage analyze results:"
log "    $PYDAMAGE_DIR"
log "  pydamage plots:"
log "    $PLOT_DIR"
log "  logs:"
log "    $LOG_DIR"
