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
| `--anchor` | ❌ | `endpoints` | Point of the mapped target distances are measured from: `endpoints` (legacy) or `center` |
| `--anchor_pos` | ❌ | `31` | 1-based position within the target used as the measurement point when `--anchor center` |
| `--nproc` | ❌ | `1` | Number of parallel worker processes (on PBS, match your `ncpus` allocation) |
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
| `location` | Where the aligned target sits: `inside` (fully within one gene), `partial` (overlaps a gene edge) or `intergenic`. Descriptive only — every class is analysed |
| `seq_start` / `seq_end` | Normalised (ascending) input coordinate window |
| `align_start` / `align_end` | Absolute genomic coordinates of the best alignment, mapped back to forward-strand genome coordinates (1-based inclusive) |
| `anchor` | Genomic coordinate the distances were measured from (1-based) |
| `anchor_source` | `aligned` (anchor mapped cleanly through the alignment), `midpoint_fallback` (anchor base gapped/uncovered, alignment midpoint used instead), or `endpoint_min` under `--anchor endpoints` |
| `anchor_in_gene` | Locus tag(s) of the gene(s) whose body covers the measurement point, `;`-separated; empty when it is intergenic. Compare with `up_gene`/`down_gene` to see when a flanking boundary belongs to a neighbouring gene |
| `which_boundary_used` | Which alignment endpoint (`start` or `end`) was closest to a flanking gene, or `center` under `--anchor center` |
| `boundary_used` | Absolute genomic coordinate of the measurement point |
| `up_gene` / `up_boundary` / `up_dist` | Locus tag, boundary coordinate, and distance (bp) of the upstream flanking gene |
| `down_gene` / `down_boundary` / `down_dist` | Locus tag, boundary coordinate, and distance (bp) of the downstream flanking gene |
| `score` | Smith–Waterman alignment score |
| `strand` | Alignment strand (`+` or `-`) |
| `error` | Error message if processing failed; `None` on success |

## **`distances.csv`**

Two-column file (`up_dist`, `down_dist`) for successfully processed entries.\
Convenient for downstream statistical analysis. Errors and targets with no\
flanking gene annotation at all are excluded; targets inside a gene are\
**included**.

## **`distances_random.csv`**

Written whenever `--n_random > 0`. Columns `accession`, `up_dist`, `down_dist` —
one row per accepted random placement, using the same anchor as the real data.
The `accession` label allows the null to be collapsed per genome the same way
the real data can be, which is what guards against pseudoreplication when many
targets come from the same genome.

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

-   **Flanking gene logic:** Under the default `--anchor endpoints`, both the\
    `start` and `end` coordinates are compared against all gene boundaries and the\
    endpoint with the shorter minimum flanking distance is reported. Note that this\
    is a min-of-two statistic and therefore shifts the reported distribution\
    downward relative to any single fixed reference point.

-   **Targets inside genes are counted.** Every target is measured to the nearest\
    gene boundary upstream and downstream of it — each gene judged in its own\
    orientation — whether or not the target sits inside a gene, and whether that\
    boundary belongs to the gene it sits in or to a neighbour. For a target inside\
    a + strand gene, for example, the downstream stop codon is usually that gene's\
    own and the upstream stop codon is the previous gene's. Earlier versions\
    dropped every target fully inside a gene, which removed the sites furthest\
    from any codon and pulled cumulative distance distributions towards zero. The\
    random null keeps inside-gene placements for the same reason. Only targets\
    with no flanking gene in either direction (GenBank records without gene/CDS\
    features) are excluded.

-   **Origin-spanning genes:** a gene that crosses the origin of a circular\
    molecule (`join(4044092..4044757,1..936)`) is handled as its real segments,\
    with its start and stop codons taken from its first and last parts in\
    biological order. Earlier versions used the feature's min/max extent, which\
    made that one gene cover the entire chromosome — so every target in any\
    genome containing such a gene was classified inside a gene and dropped — and\
    placed its start/stop codon at the genome ends.

-   **Circular molecules:** for records whose GenBank topology is `circular`,\
    distances are measured around the origin, so a target near either end of the\
    sequence reaches the next gene across it. Linear records are not wrapped.

-   **`--anchor center`:** measures instead from one fixed position inside the\
    target (`--anchor_pos`, default 31 — the first base past the midpoint of a\
    60 bp target), which for a target window centred on an insertion site is the\
    insertion point itself. The position is interpreted in the target's own\
    stranded orientation and mapped through the alignment, so gaps and\
    reverse-strand entries land on the intended nucleotide. When `--n_random > 0`\
    the random-placement null uses the same anchor offset, so real and null remain\
    directly comparable. Rows where the anchor base is not covered by the local\
    alignment fall back to the alignment midpoint and are flagged in\
    `anchor_source`.

-   **Coordinate convention:** `blastdbcmd -range` is 1-based inclusive while\
    BioPython feature coordinates are 0-based, so all distance arithmetic is done\
    0-based internally and coordinates are converted to 1-based inclusive only on\
    output. A distance of 0 therefore means the anchor nucleotide *is* the gene\
    boundary nucleotide. (Runs produced before this fix carried a +1 bp offset on\
    `up_dist` and a −1 bp offset on `down_dist`; do not compare across the two at\
    single-bp resolution.)

-   **FASTA headers:** the accession/coordinate split takes the **last** `:` or\
    `_` in the header, so accessions that themselves contain an underscore\
    (`NC_055040`, `NM_001126745`, …) are parsed rather than dropped.

-   **Error handling:** Failures (missing GenBank file, no alignment found, etc.)\
    are caught per-region and written to the output CSV with the `error` column\
    populated, so the pipeline does not abort on partial failures.

-   **Temp files:** Intermediate `blastdbcmd` FASTA extracts are written to\
    `<output_dir>/blast_tmp/` and deleted immediately after parsing. Ensure this\
    path is on a filesystem that handles concurrent writes if using `--nproc > 1`.
