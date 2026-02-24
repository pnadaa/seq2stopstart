# seq2startstop

A pipeline for mapping genomic sequences to a BLAST database, aligning them
against their reference regions, and identifying the nearest flanking gene
boundaries (start or stop codons) on both sides of each alignment.

Designed for high-throughput analysis of insertion sequences or other mobile
genetic elements across many bacterial genomes.

---

## Overview

For each sequence in an input FASTA file, this tool:

1. Parses the genomic coordinates from the FASTA header (`ACCESSION_START-END`)
2. Extracts the reference window from a local BLAST database using `blastdbcmd`
3. Performs strand-aware Smith–Waterman local alignment of the query against
   the reference window using `Bio.Align.PairwiseAligner`
4. Parses the corresponding GenBank file to retrieve gene feature boundaries
5. Identifies the closest upstream and downstream gene (by stop codon or start
   codon, configurable) relative to the alignment endpoints
6. Writes results to CSV and optionally plots a distance distribution histogram

---

## Requirements

### Software

- Python ≥ 3.10 (uses `X | Y` union type hints)
- NCBI BLAST+ (`blastdbcmd` must be on your `PATH`)

### Python Dependencies
biopython\
matplotlib


## Install with:

```bash
pip install biopython matplotlib ```

------------------------------------------------------------------------

## **Input Files**

## **FASTA file (`--fasta`)**

Each record must have a header in the format:

```text  
`>ACCESSION_START-END`

-   `ACCESSION` — the sequence accession as it appears in your BLAST database

-   `START` and `END` — 1-based genomic coordinates of the region

-   If `START > END`, the entry is treated as reverse-strand; coordinates are\
    automatically normalised to ascending order before querying `blastdbcmd`
```
**Example headers:**

```text
`>CP000253_1234000-1235500
>AE006468_500000-498200`
```
## **BLAST database (`--blast_db`)**

A local BLAST nucleotide database built with `makeblastdb`. Pass the database\
**prefix** (i.e., the value given to `makeblastdb -out`), not the directory\
containing it.

**Correct:** `--blast_db /data/blastdb/staph_genomes`\
**Incorrect:** `--blast_db /data/blastdb/`

## **GenBank files (`--genbank_dir`)**

One `.gbff` file per genome. By default, files are expected to be named\
`{ACCESSION}.gbff`. This can be customised with `--gb_naming`.

------------------------------------------------------------------------

## **Usage**

```bash
`python seq2startstop.py \
    --fasta sequences.fasta \
    --genbank_dir /path/to/gbff/ \
    --blast_db /path/to/blastdb/mydb \
    --output_dir results/`
```
## **Full argument reference**

| **Argument** | **Required** | **Default** | **Description** |
|:---|:---|:---|:---|
| `--fasta` | ✅ | — | Input FASTA file with `ACCESSION_START-END` headers |
| `--genbank_dir` | ✅ | — | Directory containing GenBank (`.gbff`) files |
| `--blast_db` | ✅ | — | BLAST+ database prefix (e.g. `/data/blastdb/mydb`), not a directory |
| `--output_dir` | ✅ | — | Directory for all output files (created if absent) |
| `--gb_naming` | ❌ | `accession` | GenBank filename pattern; use `{accession}` as placeholder, e.g. `{accession}_genomic.gbff` |
| `--boundary_type` | ❌ | `stop` | Gene boundary to measure distance to: `stop` or `start` |
| `--nproc` | ❌ | `1` | Number of parallel worker processes (on PBS, match your `ncpus` allocation) |
| `--csv_out` | ❌ | `coordinates_with_genes.csv` | Filename for the main results CSV |
| `--dist_out` | ❌ | `distances.csv` | Filename for the distances-only CSV |
| `--plot` | ❌ | off | If set, generate a histogram of upstream/downstream distances |
| `--plot_name` | ❌ | `distance_distribution.png` | Filename for the histogram image |
| `--xlim` | ❌ | off | Set the limits for the x axis - usage: --xlim 0 50 for 0-50 |
| `--verbose` | ❌ | off | Print `blastdbcmd` commands and additional progress info |

------------------------------------------------------------------------

## **Output Files**

## **`coordinates_with_genes.csv` (main output)**

One row per input sequence. Columns:

| **Column** | **Description** |
|:---|:---|
| `accession` | Genome accession |
| `query` | Original FASTA header |
| `is_reverse` | Whether the input coordinates implied a reverse-strand entry |
| `seq_start` / `seq_end` | Normalised (ascending) input coordinate window |
| `align_start` / `align_end` | Absolute genomic coordinates of the best alignment, mapped back to forward-strand genome coordinates |
| `which_boundary_used` | Which alignment endpoint (`start` or `end`) was closest to a flanking gene |
| `boundary_used` | Absolute genomic coordinate of the selected alignment endpoint |
| `up_gene` / `up_boundary` / `up_dist` | Locus tag, boundary coordinate, and distance (bp) of the upstream flanking gene |
| `down_gene` / `down_boundary` / `down_dist` | Locus tag, boundary coordinate, and distance (bp) of the downstream flanking gene |
| `score` | Smith–Waterman alignment score |
| `strand` | Alignment strand (`+` or `-`) |
| `error` | Error message if processing failed; `None` on success |

## **`distances.csv`**

Two-column file (`up_dist`, `down_dist`) for successfully processed entries.\
Convenient for downstream statistical analysis (errors are excluded).

## **`distance_distribution.png` (optional)**

Overlapping histogram of upstream and downstream distances to the configured\
boundary type, generated when `--plot` is passed.

------------------------------------------------------------------------

## **Examples**

## **Basic run (stop codon distances)**

```bash
python seq2startstop.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db /data/blastdb/staph_genomes \
    --output_dir out/`
```
## **Start codon distances, parallelised, with plot**

```bash
python seq2startstop.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db /data/blastdb/staph_genomes \
    --output_dir out/ \
    --boundary_type start \
    --nproc 8 \
    --plot \
    --xlim 0 50 \
    --plot_name start_codon_distances.png`
```
## **Custom GenBank filename pattern**

```bash
python seq2startstop.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db /data/blastdb/staph_genomes \
    --output_dir out/ \
    --gb_naming "{accession}_genomic.gbff"`
```
## **PBS job script example**

```bash
#PBS -l ncpus=16,mem=32gb,walltime=4:00:00
cd $PBS_O_WORKDIR
python seq2startstop.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db /data/blastdb/staph_genomes \
    --output_dir out/ \
    --nproc 16 \
    --plot`
```
------------------------------------------------------------------------

## **Implementation Notes**

-   **Alignment engine:** Uses `Bio.Align.PairwiseAligner` in local (Smith–Waterman)\
    mode with scoring: match `+2`, mismatch `−1`, gap open `−2`, gap extend `−0.5`.\
    Both the query and its reverse complement are aligned against the forward\
    reference; the higher-scoring strand is reported.

-   **Strand handling:** Reverse-strand entries (where `START > END` in the FASTA\
    header) are detected and coordinate-normalised before querying `blastdbcmd`.\
    Alignment coordinates are always reported in forward-strand genome space.

-   **Flanking gene logic:** For each alignment, both the `start` and `end`\
    coordinates are compared against all gene boundaries; the endpoint with the\
    shorter minimum flanking distance is used for reporting.

-   **Error handling:** Failures (missing GenBank file, no alignment found, etc.)\
    are caught per-region and written to the output CSV with the `error` column\
    populated, so the pipeline does not abort on partial failures.

-   **Temp files:** Intermediate `blastdbcmd` FASTA extracts are written to\
    `<output_dir>/blast_tmp/` and deleted immediately after parsing. Ensure this\
    path is on a filesystem that handles concurrent writes if using `--nproc > 1`.
