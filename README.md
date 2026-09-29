# SEDNA

SEDNA is an integrated command-line workflow for sedimentary ancient DNA (sedaDNA)
and environmental DNA analysis. Its core capabilities include:

1. Reference database construction (NCBI nuclear genomes + organelle sequences)
2. Bowtie2 alignment
3. ngsLCA taxonomic annotation
4. Per-taxon damage analysis (bamdam)
5. Damage curve visualization across taxonomic ranks
6. Read-count summaries per taxonomic rank
7. Per-taxon BAM extraction with coverage computation
8. Comprehensive sedaDNA authenticity evaluation
9. *de novo* assembly + BWA mapping + pydamage analysis
10. Extraction of ancient-like reads from damaged contigs
11. Benchmark wrapper for runtime / memory / tool-version profiling

---

## 0. Package Layout and Directory Conventions

The delivered package contains two sibling directories:

```
sedna_package/
├── sedna/
│   ├── bin/                  CLI entry points (sedna, sedna-*)
│   ├── scripts/              implementation scripts called by bin/
│   │                           + run_test.sh — end-to-end smoke test (section 12)
│   ├── supplementary_scripts/not called by the CLI; unsupported/optional helpers
│   ├── envs/                 sedna_unified_env.yml
│   ├── examples/             acc_list.txt, species.txt,
│   │                           ERR10878159_clean_downsample_10pct.fq.gz
│   └── README.md
└── taxdmp/                   taxonomy database — download yourself (section 3)
```

Throughout this document the examples use the following placeholder directories.
**Replace them with absolute paths on your system** (recommended on HPC, where the
working directory of compute nodes may differ from the login node):

| Placeholder | Meaning |
|---|---|
| `/project` | Root of your working area |
| `/project/ref_cave` | Reference database directory, produced by `sedna build-db` |
| `/project/taxdmp` | Taxonomy dump directory (section 3) |
| `/project/results/` | Output directory |

> **Note:** the reference database (`ref_cave`) is **not** shipped with the package.
> You must build it once with `sedna build-db` before running `sedna run`.
> The taxonomy database is also **not** shipped — download it as described in section 3.

---

## 1. Installation

### 1.1 Create the environment

The environment file `envs/sedna_unified_env.yml` is exported for **Linux x86_64**
(it pins `ld_impl_linux-64` and other Linux-only builds), so create it on the HPC,
not on macOS/Windows:

```bash
conda env create -f envs/sedna_unified_env.yml
conda activate sedna
```

### 1.2 Make the CLI executable

The scripts in `bin/` are delivered without the executable bit. Grant it once:

```bash
chmod +x /project/sedna_package/sedna/bin/sedna /project/sedna_package/sedna/bin/sedna-*
```

### 1.3 Add the CLI to `PATH`

```bash
export PATH="/project/sedna_package/sedna/bin:$PATH"
sedna --version
```

### 1.4 Verify the dependencies

```bash
which bowtie2 samtools ngsLCA bamdam python3 datasets seqkit bwa megahit pydamage
```

Additional requirements per stage:

| Stage | Requires |
|---|---|
| `build-db` | **network access**, `curl`, `datasets` (or `ncbi-genome-download`), `seqkit`, `bowtie2-build` |
| `run` | `bowtie2`, `samtools`, `ngsLCA`, `python3` |
| `damage` | `bamdam`, `samtools`, `python3` (with `pandas`) |
| `assemble-pydamage` | `megahit`, `bwa`, `samtools`, `pydamage` |

`curl` is only used to download taxonomy files, and is normally provided by the
system; the conda environment itself ships `libcurl` (library) only.

---

## 2. Prepare Input Files

Ready-to-use examples are in `sedna/examples/`:

- `acc_list.txt` — nuclear genome accessions
- `species.txt` — species/taxid list for organelle retrieval
- `ERR10878159_clean_downsample_10pct.fq.gz` — test reads (single-end, gzipped FASTQ)

### 2.1 Accessions file (nuclear genomes)

One GCA/GCF accession per line:

```
GCA_000001405.29
GCF_000001635.27
```

### 2.2 Species file (for organelle genome download)

Tab- or space-delimited, one `species + taxid` per line:

```
Homo sapiens    9606
Mus musculus    10090
```

---

## 3. Taxonomy Database (taxdmp) — Download It Yourself

The taxonomy database is **not included** in the package. Download it from the NCBI
taxonomy FTP site:

**`https://ftp.ncbi.nih.gov/pub/taxonomy/`**

SEDNA needs exactly these three files in **one directory** (e.g. `/project/taxdmp`):

| File | Source |
|---|---|
| `names.dmp` | extracted from `taxdump.tar.gz` |
| `nodes.dmp` | extracted from `taxdump.tar.gz` |
| `nucl_gb.accession2taxid` | decompressed from `nucl_gb.accession2taxid.gz` |

```bash
# Pick a permanent location for the taxonomy database
TAXDMP_DIR=/project/taxdmp
mkdir -p "$TAXDMP_DIR" && cd "$TAXDMP_DIR"

# 1) names.dmp + nodes.dmp are shipped together in taxdump.tar.gz
#    (use `curl -O` if curl is available; otherwise `wget`)
curl -O https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz names.dmp nodes.dmp
rm -f taxdump.tar.gz

# 2) accession-to-taxid mapping: the pipeline only needs nucl_gb (GenBank).
#    nucl_wgs (WGS contigs) is LARGE (100+ GB) and NOT required by SEDNA — do not
#    download it unless you have a specific reason to use it outside the pipeline.
curl -O https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
# wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip -c nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid
rm -f nucl_gb.accession2taxid.gz
```

> **Important**
> - The three file names above are **hard-coded** in the pipeline: `names.dmp`,
>   `nodes.dmp`, `nucl_gb.accession2taxid`. Do not rename them.
> - `nucl_gb.accession2taxid` must be the **uncompressed** file; a `.gz` file will
>   not be recognised.
> - The full `nucl_gb.accession2taxid` is large (~20+ GB uncompressed). Make sure the
>   target filesystem has enough space.

---

## 4. Reference Database Construction

```bash
sedna build-db \
  -a /project/sedna_package/sedna/examples/acc_list.txt \
  -o /project/ref_cave \
  --species-file /project/sedna_package/sedna/examples/species.txt \
  --threads 16
```

Output:

- `/project/ref_cave/ncbi_reference.*.bt2` (or `*.bt2l`) — Bowtie2 index
- `/project/ref_cave/work_build_db/combined.dedup.fasta` — merged reference FASTA

Options:

| Argument | Description | Default |
|---|---|---|
| `-a, --accessions FILE` | NCBI assembly accessions (`GCA_`/`GCF_`), one per line | — |
| `-o, --outdir DIR` | Output directory (**required**) | — |
| `--species-file FILE` | Two-column species/taxid list (organelle retrieval) | `./species.txt` |
| `--use-existing-genomes-only` | Skip download; merge local FASTAs from `--genome-dir` | off |
| `--genome-dir DIR` | Local nuclear FASTA directory (recursive: `.fa/.fasta/.fna`, optionally `.gz`) | — |
| `--api-key KEY` / `--email EMAIL` | NCBI credentials (raise rate limits) | — |
| `--chunk-size N` | Download chunk size | `50` |
| `--threads N` | Threads for dedup / `bowtie2-build` | `8` |
| `--force-stage S` | Re-run from a stage: `genome\|organelle\|combine\|dedup\|collapse\|bowtie2\|all` | — |
| `--collapse-by-taxon` | Concatenate scaffolds of the same assembly into one record (joined by 100 Ns) | off |
| `--index-prefix PREFIX` | Bowtie2 index prefix | `OUTDIR/ncbi_reference` |

---

## 5. Main Workflow: Alignment + Classification

```bash
sedna run \
  -i /project/sedna_package/sedna/examples/ERR10878159_clean_downsample_10pct.fq.gz \
  -o /project/results/sample_test \
  --ref-dir /project/ref_cave \
  --taxdb-dir /project/taxdmp \
  -t 16
```

Main outputs (in `-o` directory):

- `out.name_sorted.bam`
- `lca_results.lca` (or another `lca_results*` file)
- `bowtie2.log`, `ngsLCA.stdout.log`, `ngsLCA.stderr.log`

Options:

| Argument | Description | Default |
|---|---|---|
| `-i, --input FILE` | Input FASTQ / FASTQ.gz (**required**) | — |
| `-o, --output DIR` | Output directory (**required**) | — |
| `--ref-dir DIR` | Reference directory produced by `sedna build-db` (**required**) | — |
| `--taxdb-dir DIR` | Directory with `names.dmp`, `nodes.dmp`, `nucl_gb.accession2taxid` (**required**) | — |
| `-t, --threads INT` | Threads | `8` |
| `--skip-bowtie2` | Reuse an existing `out.name_sorted.bam` | off |
| `--skip-ngslca` | Reuse an existing `lca_results*` | off |
| `--skip-abundance` | Deprecated no-op (abundance output is disabled) | — |

> **Note:** This version does **not** generate abundance tables or abundance plots.

---

## 6. Damage Analysis (bamdam)

`sedna damage` runs bamdam compute on one or more taxonomic ranks, extracts all taxa
from the bamdam results, and optionally performs BAM extraction, coverage calculation,
and comprehensive evaluation.

### 6.1 Multi-rank run (recommended)

```bash
sedna damage \
  --sedna-dir /project/results/sample_test \
  --ranks "species,genus,family" \
  --stranded ds \
  --outdir /project/results/sample_test_bamdam
```

For each rank (`species`, `genus`, `family`), the pipeline runs these steps:

| Step | Description |
|------|-------------|
| [1/5] | bamdam compute |
| [2/5] | Generate read-count summary table (`read_counts.<rank>.tsv`) |
| [3/5] | Generate damage plots (`bamdam plotdamage` smiley plots) |
| [4/5] | Extract all taxa from the bamdam TSV |
| [5/5] | Extract BAM + coverage (only when `--extract-bam` is set) |

### 6.2 Output structure

```
<outdir>/
├── species/
│   ├── compute/
│   │   ├── bamdam.species.upto_species.tsv
│   │   └── bamdam.species.upto_species.subs
│   ├── read_counts.species.tsv
│   ├── damage_plots/
│   │   ├── 9606_Homo_sapiens_damage.png
│   │   └── ...
│   ├── reads_species/          (with --extract-bam)
│   ├── coverage_analysis/      (with --extract-bam)
│   └── logs/
├── genus/
│   └── ...
├── family/
│   └── ...
├── comprehensive_adna_report.tsv      (with --evaluate)
└── comprehensive_adna_detailed.tsv    (with --evaluate)
```

### 6.3 Read-count summary table

`read_counts.<rank>.tsv` is generated automatically for every rank. Columns:

| Column | Description |
|--------|-------------|
| `TaxNodeID` | NCBI taxonomy ID |
| `TaxName` | Scientific name |
| `Rank` | Taxonomic rank |
| `TotalReads` | Total reads assigned to this taxon |
| `UnaggregatedReads` | Reads before aggregation |

Sorted by `TotalReads` descending.

### 6.4 Damage plots

Per-taxon "smiley" plots are generated using `bamdam plotdamage` from the `.subs` file.
Each plot shows the canonical postmortem damage pattern (5' C→T and 3' G→A substitution
frequencies by position) for a single taxon.

Plots are saved as PNG files in `<outdir>/<rank>/damage_plots/`.

---

## 7. Per-taxon BAM Extraction and Coverage

Add `--extract-bam` to extract per-taxon BAM files and compute coverage:

```bash
sedna damage \
  --sedna-dir /project/results/sample_test \
  --ranks "species,genus,family" \
  --stranded ds \
  --extract-bam \
  --outdir /project/results/sample_test_bamdam
```

BAM output and coverage files:

- Per-rank BAM: `<outdir>/<rank>/reads_<rank>/<taxid>.bam`
- Coverage: `<outdir>/<rank>/coverage_analysis/<taxid>_coverage.txt`
- Coverage summary: `<outdir>/<rank>/coverage.<rank>.tsv`

---

## 8. Comprehensive sedaDNA Authenticity Evaluation

After running bamdam analysis with `species`, `genus`, and `family` ranks, add
`--evaluate` to perform a comprehensive evaluation of aDNA authenticity. This combines
three dimensions:

1. **Terminal damage pattern** (35 points): 5' C→T damage, decay pattern, 5'/3' asymmetry
2. **Coverage uniformity + read support** (30 points): coverage breadth, Gini coefficient, read count
3. **Taxonomic specificity** (35 points): read-sink detection, LCA degradation, ANI, GC bias

```bash
sedna damage \
  --sedna-dir /project/results/sample_test \
  --ranks "species,genus,family" \
  --stranded ds \
  --extract-bam \
  --evaluate \
  --outdir /project/results/sample_test_bamdam
```

> `--evaluate` requires the `species`, `genus`, and `family` TSVs to all be present.
> If any rank is missing, the evaluation step is skipped.

Output files:

- `comprehensive_adna_report.tsv` — simplified report, top species sorted by score
- `comprehensive_adna_detailed.tsv` — detailed report with all species and intermediate metrics

Confidence levels:

- >= 65: **HIGH** — High-confidence ancient DNA
- >= 45: **MEDIUM** — Medium confidence
- >= 30: **LOW** — Low confidence, needs further validation
- >= 12: **VERY_LOW** — Very low confidence, likely false positive
- < 12: **BACKGROUND** — Background noise or contamination

---

## 9. Optional Arguments for `sedna damage`

| Argument | Description | Default |
|----------|-------------|---------|
| `--sedna-dir DIR` | `sedna run` output directory (**required**); must contain `out.name_sorted.bam` and an `lca_results*` file | — |
| `--ranks` | Comma-separated taxonomic ranks | `species` |
| `--stranded` | Library strandedness: `ss` or `ds` | `ds` |
| `--upto` | Aggregation rank for bamdam | Same as rank |
| `--mode` | bamdam compute mode | `1` |
| `--udg` | Enable UDG-treated library mode | off |
| `--ref-fasta` | Reference FASTA for calmd (memory-intensive) | — |
| `--k` | k-mer length | `29` |
| `--safe-bam` | Apply safe BAM filtering | off |
| `--extract-bam` | Extract per-taxon BAM and compute coverage | off |
| `--evaluate` | Run comprehensive authenticity evaluation | off |
| `--outdir` | Output directory | `bamdam_native` |

---

## 10. Extended Subcommands

Besides `build-db` / `run` / `damage`, the CLI provides three more subcommands.
All arguments below are **positional** (not `--flags`), and output is written relative
to the current working directory, so run them from a writable directory.

### 10.1 `assemble-pydamage` — assembly + mapping + pydamage

```bash
sedna assemble-pydamage <R1.fastq[.gz]> <R2.fastq[.gz]|-> <sample_name> <threads>
```

Single-end example (use `-` as the R2 placeholder):

```bash
sedna assemble-pydamage XD_S_136_lt100bp.fastq - sample136 16
```

Paired-end example:

```bash
sedna assemble-pydamage sample_R1.fastq.gz sample_R2.fastq.gz sample136 16
```

Steps: MEGAHIT assembly → BWA-MEM mapping to contigs → `samtools calmd` (MD tag) →
`pydamage analyze`.

Outputs under `<sample_name>.work/`:

| Path | Content |
|---|---|
| `assembly/final.contigs.fa` | MEGAHIT contigs |
| `mapping/<sample>.sorted.bam` | sorted alignments |
| `mapping/<sample>.sorted.MD.bam` | BAM with MD tag (pydamage input) |
| `pydamage/pydamage_results/pydamage_results.csv` | pydamage results |

> **Warning:** if `<sample_name>.work` already exists it is **deleted and rebuilt**,
> because MEGAHIT requires its `-o` target not to pre-exist.

### 10.2 `extract-adna` — extract ancient-like reads

```bash
sedna extract-adna <pydamage_results.csv> <sorted.bam> <out_prefix> <single|paired> [threads]
```

Example:

```bash
sedna extract-adna \
  sample136.work/pydamage/pydamage_results/pydamage_results.csv \
  sample136.work/mapping/sample136.sorted.MD.bam \
  sample136.ancient_like single 16
```

High-confidence contigs are selected with fixed thresholds:
`predicted_accuracy >= 0.95`, `qvalue <= 0.05`, `nb_reads_aligned >= 100`,
`coverage >= 3`, `CtoT-0 >= 0.03` (edit the script to change them). Reads aligned to
those contigs are then extracted and converted back to FASTQ.

Output directory `<out_prefix>.extract_ancient_like/` contains the contig list, BED,
filtered BAM, and FASTQ (single-end: `<prefix>.ancient_like.fastq.gz`; paired-end:
`.R1/.R2/.singleton/.other.fastq.gz`).

### 10.3 `bench` — record runtime / memory / versions

```bash
# Wrap a command (note the required "--" separator)
sedna bench --label run -- sedna run -i reads.fq.gz -o out --ref-dir ref --taxdb-dir tax -t 16

# Aggregate all records into a Markdown table
sedna bench --summary ./bench
```

| Argument | Description | Default |
|---|---|---|
| `--label LABEL` | Output file label (default: inferred from the command) | — |
| `--out-dir DIR` | Directory for JSON records | `./bench` |
| `--summary DIR` | Summary mode: print a Markdown table for all `*.json` in `DIR` | — |

---

## 11. `bin/` vs `scripts/` vs `supplementary_scripts/`

- **`bin/`** contains the official CLI entry points. Always call SEDNA through
  `sedna <subcommand>`; do not edit or invoke `bin/sedna-*` directly.
- **`scripts/`** contains the implementations and the end-to-end smoke test:
  - `run_sedna.sh` — called by `sedna run`
  - `evaluate_adna_comprehensive.py` — called by `sedna damage --evaluate`
  - `plot_damage_curves.py` — called by `sedna damage` (damage plots)
  - `generate_read_counts.py` — called by `sedna damage` (read-count tables)
  - `run_test.sh` — standalone smoke test (section 12)
  - `filter_bamdam_results.py` — optional post-processing helper, not called by the CLI
- **`supplementary_scripts/`** contains optional, unsupported helper scripts that
  are not wired into the CLI. They are preserved for reproducibility of ad-hoc
  analyses, but are not part of the supported workflow.

---

## 12. End-to-End Validation (Smoke Test)

The package ships with a self-contained validation script that exercises the three
core stages on the example data in `examples/`.  It uses:

- `examples/acc_list.txt`
- `examples/species.txt`
- `examples/ERR10878159_clean_downsample_10pct.fq.gz`

### 12.1 Quick-start: install env + launch the example

Run these commands **from inside the `sedna/` directory** of the unpacked package
on a Linux HPC node with internet access (required to download the NCBI taxdmp
on first use):

```bash
# 1) Install the pinned conda environment
conda env create -f envs/sedna_unified_env.yml
conda activate sedna

# 2) Run the end-to-end example in the background with nohup
#    (stderr and stdout both go to test_ERR10878159_10pct.log)
nohup bash scripts/run_test.sh \
    > test_ERR10878159_10pct.log  2>&1  &
```

You can follow progress with:

```bash
tail -f test_ERR10878159_10pct.log
```

### 12.2 How to confirm the test finished successfully

Open or `tail -n 40 test_ERR10878159_10pct.log`.  When the test completes
without errors, the last block of output will contain exactly the banner:

```
=======================================================
SEDNA smoke test PASSED  (YYYY-MM-DD HH:MM:SS)
=======================================================
```

for example:

```
=======================================================
SEDNA smoke test PASSED  (2026-09-29 11:39:13)
=======================================================
```

If you see **`SEDNA smoke test PASSED`**, the full pipeline (download
taxdmp → `sedna build-db` → `sedna run` → `sedna damage`) has completed on the
shipped example data and every required output file was produced and non-empty.

Immediately after the banner the log also lists the absolute paths to all key
artifacts (merged reference FASTA, BAM, LCA results, bamdam tables, damage
plots) and the `rm -rf` command needed to clean up the temporary test work
directory.

### 12.3 Reusing an existing taxonomy database

If you already have a local copy of the NCBI taxonomy dump (`names.dmp`,
`nodes.dmp`, and the **uncompressed** `nucl_gb.accession2taxid` all in one
directory), pass its absolute path as the first positional argument to avoid
re-downloading it every run:

```bash
conda activate sedna
cd /project/sedna_package/sedna

# Reuse /project/taxdmp (fast, no network for taxdmp; still needs network for
# sedna build-db to fetch the reference sequences).
nohup bash scripts/run_test.sh /project/taxdmp \
    > test_ERR10878159_10pct.log  2>&1  &
```

### 12.4 What `run_test.sh` does step-by-step

1. Grants the executable bit to every entry in `bin/` and adds it to `PATH`.
2. Verifies `sedna`, `bowtie2`, `samtools`, `ngsLCA`, `bamdam`, `python3` are on `PATH`.
3. Downloads or reuses the NCBI taxonomy dump under a unique timestamped work dir
   (`sedna/test_YYYYMMDD_HHMMSS/`).
4. **Stage 1** `sedna build-db` → asserts merged FASTA + Bowtie2 index exist.
5. **Stage 2** `sedna run` → asserts `out.name_sorted.bam` and an `lca_results*` file exist.
6. **Stage 3** `sedna damage` (species rank, `--stranded ds`) → asserts bamdam
   TSV/subs, `read_counts.species.tsv`, and at least one damage plot PNG exist.
7. On success prints the paths of all produced artifacts and how to clean up.

All output is written to the unique test work directory so successive runs never
clobber each other; a combined log is saved to `test_<stamp>/test.log`.

---

## 13. One-click Test Workflow

```bash
PKG=/project/sedna_package

# 1) Build the reference database (needs network)
sedna build-db \
  -a $PKG/sedna/examples/acc_list.txt \
  -o /project/ref_cave \
  --species-file $PKG/sedna/examples/species.txt \
  --threads 16

# 2) Main workflow
sedna run \
  -i $PKG/sedna/examples/ERR10878159_clean_downsample_10pct.fq.gz \
  -o /project/results/sample_test \
  --ref-dir /project/ref_cave \
  --taxdb-dir /project/taxdmp \
  -t 16

# 3) Damage analysis + read counts + damage plots + BAM extraction + evaluation
sedna damage \
  --sedna-dir /project/results/sample_test \
  --ranks "species,genus,family" \
  --stranded ds \
  --extract-bam \
  --evaluate \
  --outdir /project/results/sample_test_bamdam
```

---

## 14. FAQ

**Q: `sedna: permission denied`.**
The `bin/` scripts need the executable bit; run the `chmod` command in section 1.2.

**Q: `sedna run` stops with `names.dmp does not exist or is empty`.**
The taxonomy database is not shipped. Download it as described in section 3 and point
`--taxdb-dir` to the directory containing `names.dmp`, `nodes.dmp`, and the
**uncompressed** `nucl_gb.accession2taxid`.

**Q: `Could not find bowtie2 index`.**
The reference database has not been built yet. Run `sedna build-db` first (section 4),
and pass its `-o` directory to `--ref-dir`.

**Q: No `lca_results*` found.**
Make sure `sedna run` completed successfully and generated `lca_results*` in the output directory.

**Q: Why are there no abundance tables or plots?**
Abundance output is intentionally disabled in this version. `sedna damage` and the
comprehensive evaluation use bamdam outputs directly.

**Q: `bamdam` fails with `StopIteration`.**
The pipeline automatically cleans LCA read-ID suffixes and filters to primary alignments.
If the issue persists, check whether the input BAM has been externally modified.

**Q: Can I filter taxa by read count or select top N?**
The current version processes all taxa from bamdam results without read-count or top-N
filtering. Use the generated `read_counts.<rank>.tsv` to inspect and filter taxa downstream.
