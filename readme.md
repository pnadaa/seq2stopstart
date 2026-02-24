# IS Boundary Locator

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
   the reference window
4. Parses the corresponding GenBank file to retrieve gene feature boundaries
5. Identifies the closest upstream and downstream gene (by stop codon or start
   codon, configurable) relative to the alignment
6. Writes results to CSV and optionally plots a distance distribution histogram

---

## Requirements

### Software

- Python ≥ 3.10 (uses `X | Y` union type hints)
- NCBI BLAST+ (`blastdbcmd` must be on your `PATH`)

### Python Dependencies
biopython\
matplotlib


Install with:

```bash
pip install biopython matplotlib ```

> **⚠️ Deprecation Warning: This script uses `Bio.pairwise2`, which was\
> deprecated in Biopython 1.80. It will still function in current releases but\
> may be removed in a future version. Consider migrating to\
> `Bio.Align.PairwiseAligner` for long-term compatibility.**

------------------------------------------------------------------------

## **Input Files**

## **FASTA file (`--fasta`)**

Each record must have a header in the format:

```text

`>ACCESSION_START-END`

-   `ACCESSION` — the sequence accession as it appears in your BLAST database

-   `START` and `END` — 1-based genomic coordinates of the region

-   If `START > END`, the entry is treated as reverse-strand and coordinates are\
    automatically normalised
```

**Example headers:**

```         


`>CP000253_1234000-1235500
>AE006468_500000-498200`
```
## **BLAST database (`--blast_db_dir`)**

A local BLAST nucleotide database built with `makeblastdb`, containing the\
genomes corresponding to the accessions in your FASTA file.

## **GenBank files (`--genbank_dir`)**

One `.gbff` file per genome. By default, files are expected to be named\
`{ACCESSION}.gbff`. This can be customised with `--gb_naming`.

------------------------------------------------------------------------

## **Usage**

```         
bash
```

`python is_boundary_locator.py \
    --fasta sequences.fasta \
    --genbank_dir /path/to/gbff/ \
    --blast_db_dir /path/to/blastdb/mydb \
    --output_dir results/`

## **Full argument reference**

| **Argument** | **Required** | **Default** | **Description** |
|:---|:---|:---|:---|
| `--fasta` | ✅ | — | Input FASTA file with `ACCESSION_START-END` headers |
| `--genbank_dir` | ✅ | — | Directory containing GenBank (`.gbff`) files |
| `--blast_db_dir` | ✅ | — | Path to the BLAST database (prefix, not directory) |
| `--output_dir` | ✅ | — | Directory for all output files |
| `--gb_naming` | ❌ | `accession` | GenBank filename pattern; use `{accession}` as placeholder, e.g. `{accession}_genomic.gbff` |
| `--boundary_type` | ❌ | `stop` | Gene boundary to measure distance to: `stop` or `start` |
| `--nproc` | ❌ | `1` | Number of parallel worker processes |
| `--csv_out` | ❌ | `coordinates_with_genes.csv` | Filename for the main results CSV |
| `--dist_out` | ❌ | `distances.csv` | Filename for the distances-only CSV |
| `--plot` | ❌ | off | If set, generate a histogram of upstream/downstream distances |
| `--plot_name` | ❌ | `distance_distribution.png` | Filename for the histogram image |
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
| `align_start` / `align_end` | Absolute genomic coordinates of the best alignment |
| `which_boundary_used` | Whether `start` or `end` of the alignment was used to find the closest gene |
| `boundary_used` | Absolute coordinate of the selected alignment boundary |
| `up_gene` / `up_boundary` / `up_dist` | Locus tag, boundary coordinate, and distance of the upstream flanking gene |
| `down_gene` / `down_boundary` / `down_dist` | Locus tag, boundary coordinate, and distance of the downstream flanking gene |
| `score` | Smith–Waterman alignment score |
| `strand` | Alignment strand (`+` or `-`) |
| `error` | Error message if processing failed; `None` on success |

## **`distances.csv`**

Two-column file (`up_dist`, `down_dist`) containing only the numeric distances\
for successfully processed entries. Convenient for downstream statistical\
analysis.

## **`distance_distribution.png` (optional)**

Histogram of upstream and downstream distances to the configured boundary type\
(stop or start codons), generated when `--plot` is passed.

------------------------------------------------------------------------

## **Examples**

## **Basic run**

```         
bash
```

`python is_boundary_locator.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db_dir /data/blastdb/staph_genomes \
    --output_dir out/`

## **Measure distance to start codons, with parallelisation and a plot**

```         
bash
```

`python is_boundary_locator.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db_dir /data/blastdb/staph_genomes \
    --output_dir out/ \
    --boundary_type start \
    --nproc 8 \
    --plot \
    --plot_name start_codon_distances.png`

## **Custom GenBank filename pattern**

```         
bash
```

`python is_boundary_locator.py \
    --fasta my_IS_sequences.fasta \
    --genbank_dir /data/gbff/ \
    --blast_db_dir /data/blastdb/staph_genomes \
    --output_dir out/ \
    --gb_naming "{accession}_genomic.gbff"`

------------------------------------------------------------------------

## **Notes**

-   The temporary directory `<output_dir>/blast_tmp/` is created automatically\
    for intermediate `blastdbcmd` FASTA files and cleaned up per-query.

-   Sequences that fail (missing GenBank file, no alignment found, etc.) are\
    written to the output CSV with their `error` column populated rather than\
    halting the pipeline.

-   Alignment uses `localms` scoring: match `+2`, mismatch `−1`, gap open `−2`,\
    gap extend `−0.5`.

-   When using `--nproc > 1`, ensure the temp directory is on a filesystem that\
    handles concurrent writes safely (e.g., avoid NFS with high contention).

```         