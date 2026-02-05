# Genomic Region-Gene Distance Analysis Script

## Overview

This script analyzes a set of genomic sequences described in a FASTA file (with headers giving the genome accession and genomic coordinates). For each region, it:

-   Maps the region sequence to a corresponding genome using a local BLAST database.

-   Looks up the closest annotated upstream and downstream genes (from a local GenBank file).

-   Calculates the distance from the region to the stop codon of each flanking gene, taking gene strand into account.

-   Outputs these results as CSV files and (optionally) a summarized histogram plot.

All database and annotation access is performed locally. No internet/Entrez queries are issued.
Script logic and documentation were assisted by Gemini 3 Pro via Perplexity Pro.

## Requirements

-   Python 3.7 or higher

-   Biopython

-   matplotlib

-   BLAST+ suite (blastn, makeblastdb)

-   Local BLAST nucleotide databases (one per genomic accession)

-   Local GenBank files (.gbff/.gbk) for your genomic sequences

Install required Python packages using:\
pip install biopython matplotlib

## Easy dependency installation and running

First clone the git repository using\
`git clone https://github.com/pnadaa/seq2stopstart.git`

To run the program, install the uv package manager:\
`curl -LsSf https://astral.sh/uv/install.sh | sh`\
or\
`pip install uv`

Then run `uv sync` in the main directory to sync all packages and install dependencies.


## Input Preparation

1.  FASTA file:

    -   Each sequence header must use format:\
        \<ACCESSION\>\_\<START\>-\<END\>\
        Example:

        > **CP040883_1799184-1799626\
        > AGGGUUAAA...**

2.  GenBank files:

    -   Download GenBank files for each accession and place in the folder indicated by --genbank_dir.

    -   Files should be named as \<accession\>.gbff or according to your naming convention.

3.  BLAST databases:

    -   For each genome, prepare a BLAST nucleotide database with:\
        makeblastdb -in CP040883.fna -dbtype nucl -out ./blastdbs/CP040883

    -   The database for each accession should be in the folder specified by --blast_db_dir.

## Usage

Example usage:


`python your_script_name.py \
  --fasta regions.fasta \
  --genbank_dir ./genbanks \
  --blast_db_dir ./blastdbs \
  --output_dir ./results \
  --nproc 4 \
  --verbose \
  --plot`

Required arguments:\
--fasta Path to input FASTA file\
--genbank_dir Path to folder containing GenBank files\
--blast_db_dir Path to folder containing BLAST databases for each genome\
--output_dir Path to output directory (all results and temp files go here)

Optional arguments:\
--gb_naming Naming pattern for GenBank files (default: accession, e.g. CP040883.gbff)\
--nproc Number of processes to use (default: 1)\
--csv_out Name of full info output CSV (default: coordinates_with_genes.csv)\
--dist_out Name of distances CSV (default: distances.csv)\
--plot If present, saves histogram plot of distances as PNG file\
--verbose Enable verbose logging/status output

## Output

All outputs are placed in the directory specified by --output_dir:

-   \[csv_out\]: CSV with all region information, gene mapping, distances, alignment orientation, error reporting.

-   \[dist_out\]: CSV with just the distances, for plotting/comparative work.

-   distance_distribution.png: If --plot is set, contains a histogram plot of upstream and downstream distances.

-   blast_tmp/: Temporary folder for BLAST queries and intermediate files.

## Notes

-   Both GenBank and BLAST DBs must be present locally for every genome accession.

-   Mapping is done with blastn, checking both the sequence and its reverse complement.

-   Distances are in base pairs, measured from the boundary of each mapped region to the nearest stop codon (based on annotation strand).

-   Output CSV contains an "error" field for cases where mapping or annotation fails.

## Contact

For support or bug reports, contact the script author.\
If using on an HPC/cluster, or with very large datasets, consider running with --nproc and ensure /tmp or --output_dir has adequate space for BLAST files.