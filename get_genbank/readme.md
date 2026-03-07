## Overview

**fetch_genbank.py** is a Python script that automates the downloading of NCBI GenBank (.gbff) files using accession numbers from a multifasta file.\
It first attempts to use the NCBI Datasets CLI, and if that fails for any accession, it falls back to retrieving the sequence using Biopython’s Entrez.\
The script requires you to provide an email address for NCBI-compliant Entrez access.

## Features

-   Extracts unique NCBI accession numbers from FASTA headers formatted as **`>CP043804_723996-724438`**

-   Attempts download using NCBI Datasets CLI for each accession

-   Automatically falls back to Biopython Entrez efetch if the CLI fails

-   Renames each output GenBank file to **`<ACCESSION>.gbff`**

-   Parallel download support for fast batch operation

-   **Email input required** for NCBI Entrez access (passed via argparse)

## Requirements

-   **Python 3.6+**

-   **NCBI Datasets CLI** (**`datasets`** tool) installed and in your system PATH\
    Install via conda (recommended):

    `conda install -c conda-forge ncbi-datasets-cli`

-   **Biopython** for Entrez fallback:
-   **certifi** for Entrez fallback:

    `pip install biopython certifi`

-   Uses only standard Python libraries: **`argparse`**, **`subprocess`**, **`os`**, **`shutil`**, **`zipfile`**, **`tempfile`**, **`concurrent.futures`**

## Usage


`python fetch_genbank.py -i <input_fasta> -o <output_folder> -e <your_email> [-t N_THREADS]`

**Arguments:**

-   **`-i`**, **`--input`** : Input multifasta file

-   **`-o`**, **`--output`** : Output directory for downloaded and renamed GenBank files

-   **`-e`**, **`--email`** : Your email address (required by NCBI for Entrez)

-   **`-t`**, **`--threads`** : (Optional) Number of parallel download threads (default: 4)

-   **`-r`**, **`--resume`** : (Optional) Skip genomes which have already been downloaded

**Example:**


`python fetch_genbank.py -i example_multifasta.fasta -o downloaded_gbffs -e your@email.com -t 8 -r`

## FASTA Header Format

Headers should look like:


`>CP043804_723996-724438
>NZ_LN831024_12345-67890`

or

`>CP043804:723996-724438
>NZ_LN831024:12345-67890`

The script extracts **`CP043804`** and **`NZ_LN831024`** as accessions.

## How it works

1.  Reads the FASTA file and extracts all unique accession names.

2.  For each accession:

    -   Attempts download using NCBI Datasets CLI for GenBank record.

    -   If the CLI fails, fetches the record with Biopython's Entrez efetch (using your provided email).

3.  Outputs each as **`<ACCESSION>.gbff`** in the output folder.

## Notes

-   Accessions failing both download methods will be reported.

-   Output directory is created if not present.

-   The provided email is required by NCBI for responsible Entrez API use.

-   Supports both chromosome and WGS/refseq accessions.

## Troubleshooting

-   Make sure the NCBI Datasets CLI (**`datasets`**) is installed and available.

-   Ensure your FASTA headers match the expected format.

-   You must provide a valid email address for Entrez.

-   For proxy/network errors and API limits, see NCBI Datasets and Biopython documentation.