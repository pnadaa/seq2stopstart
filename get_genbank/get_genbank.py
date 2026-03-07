import os
import certifi
os.environ['SSL_CERT_FILE'] = certifi.where()

import argparse
import subprocess
import shutil
import zipfile
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez
from urllib.error import HTTPError

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Parse multifasta headers, download GenBank files via NCBI Datasets CLI, with Entrez fallback and rate limiting."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input multifasta file."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output directory for .gbff files."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=2,
        help="Number of parallel download threads (default: 2 for safety)."
    )
    parser.add_argument(
        "-e", "--email",
        required=True,
        help="Your email address for Entrez access (required by NCBI)."
    )
    parser.add_argument(
        "-k", "--api_key",
        required=False,
        help="NCBI API key to increase your Entrez rate limit (optional)."
    )
    parser.add_argument(
        "-r", "--resume",
        action="store_true",
        help="Skip accessions whose .gbff file already exists in the output directory."
    )
    return parser.parse_args()

def get_accessions_from_fasta(fasta_path):
    accessions = set()
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    header = line[1:].strip()
                    parts1 = header.rsplit('_', 1)
                    parts2 = header.rsplit(':', 1)
                    if len(parts1) == 2:
                        accessions.add(parts1[0])
                    elif len(parts2) == 2:
                        accessions.add(parts2[0])
                    else:
                        print(f"Warning: Could not parse header format '{header}'. Skipping.")
    except FileNotFoundError:
        print(f"Error: Input file '{fasta_path}' not found.")
        exit(1)
    return list(accessions)


def fetch_entrez_genbank(accession, output_dir, email, api_key=None,
                         sleep_interval=0.34, error_file=None):
    """
    Attempt to fetch a GenBank record via Entrez efetch.
    Retries on ANY error up to 5 total attempts, waiting 1/3/5/7 s between them.
    On final failure, appends the accession + reason to error_file (if provided).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        sleep_interval = 0.12  # up to 10 req/s with API key

    max_attempts = 5
    last_error = "Unknown error"

    for attempt in range(max_attempts):
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = handle.read()
            handle.close()
            out_path = os.path.join(output_dir, f"{accession}.gbff")
            with open(out_path, "w") as fh:
                fh.write(record)
            time.sleep(sleep_interval)
            return f"Success (Entrez): {accession}"

        except HTTPError as err:
            last_error = f"HTTP Error {err.code}: {err.reason}"
        except Exception as e:
            last_error = str(e)

        # Don't sleep after the final attempt
        if attempt < max_attempts - 1:
            wait = RETRY_WAITS[attempt]
            print(f"  Entrez attempt {attempt + 1}/{max_attempts} failed for {accession} "
                  f"({last_error}). Retrying in {wait}s...")
            time.sleep(wait)

    # All attempts exhausted — write to error file
    final_msg = f"Error: All {max_attempts} Entrez attempts failed for {accession}: {last_error}"
    if error_file:
        with open(error_file, "a") as ef:
            ef.write(f"{accession}\t{last_error}\n")
    return final_msg


def download_worker(accession, output_dir, email, api_key=None, resume=False, error_file=None):
    final_output_path = os.path.join(output_dir, f"{accession}.gbff")
    if resume and os.path.exists(final_output_path):
        return f"Skipped (already exists): {accession}"

    # Bypass CLI for nucleotide accessions — datasets CLI expects GCA_/GCF_
    if not (accession.startswith("GCA_") or accession.startswith("GCF_")):
        return fetch_entrez_genbank(accession, output_dir, email, api_key, error_file=error_file)

    with tempfile.TemporaryDirectory() as temp_dir:
        zip_filename = os.path.join(temp_dir, f"{accession}.zip")
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "gbff",
            "--filename", zip_filename
        ]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if os.path.exists(zip_filename):
                with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                    zip_ref.extractall(temp_dir)
                found_file = None
                for root, dirs, files in os.walk(temp_dir):
                    for file in files:
                        if file.endswith("genomic.gbff"):
                            found_file = os.path.join(root, file)
                            break
                if found_file:
                    shutil.move(found_file, final_output_path)
                    return f"Success: {accession}"
                else:
                    print(f"Warning: No .gbff file found for {accession}. Retrying with Entrez efetch...")
            else:
                print(f"Warning: Datasets CLI failed (zip not created) for {accession}. Retrying with Entrez efetch...")
        except subprocess.CalledProcessError:
            print(f"Warning: Datasets CLI failed for {accession}. Retrying with Entrez efetch...")
        except Exception as e:
            print(f"Error: Datasets CLI exception for {accession}: {str(e)}. Retrying with Entrez efetch...")

    # Fallback to Entrez with retry logic
    return fetch_entrez_genbank(accession, output_dir, email, api_key,
                                error_file=error_file)


def main():
    args = parse_arguments()
    if not shutil.which("datasets"):
        print("Error: 'datasets' CLI tool not found in PATH. Please install NCBI Datasets CLI.")
        exit(1)

    os.makedirs(args.output, exist_ok=True)
    error_file = os.path.join(args.output, "failed_accessions.tsv")

    print(f"Reading accessions from {args.input}...")
    accessions = get_accessions_from_fasta(args.input)
    if not accessions:
        print("No valid accessions found. Exiting.")
        exit(0)

    api_key = args.api_key if args.api_key else None
    print(f"Found {len(accessions)} unique accessions. Starting parallel download with {args.threads} threads...")
    print(f"Failed accessions will be logged to: {error_file}")

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_acc = {
            executor.submit(download_worker, acc, args.output, args.email,
                            api_key, args.resume, error_file): acc
            for acc in accessions
        }
        for future in as_completed(future_to_acc):
            print(future.result())

    print(f"\nProcessing complete. Files saved to: {args.output}")
    if os.path.exists(error_file):
        print(f"Some accessions failed — see: {error_file}")


if __name__ == "__main__":
    main()
