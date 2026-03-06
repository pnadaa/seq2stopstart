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
                    parts = header.rsplit('_', 1)
                    if len(parts) == 2:
                        acc = parts[0]
                        accessions.add(acc)
                    else:
                        print(f"Warning: Could not parse header format '{header}'. Skipping.")
    except FileNotFoundError:
        print(f"Error: Input file '{fasta_path}' not found.")
        exit(1)
    return list(accessions)

def fetch_entrez_genbank(accession, output_dir, email, api_key=None, max_retries=3, retry_wait=60, sleep_interval=0.34):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        sleep_interval = 0.12  # Can do up to 10 req/sec with API key

    retries = 0
    while retries < max_retries:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = handle.read()
            handle.close()
            out_path = os.path.join(output_dir, f"{accession}.gbff")
            with open(out_path, "w") as out_handle:
                out_handle.write(record)
            time.sleep(sleep_interval)  # Enforce rate limit
            return f"Success (Entrez): {accession}"
        except HTTPError as err:
            if err.code == 429:
                print(f"Rate limit hit for {accession}. Waiting {retry_wait} seconds...")
                time.sleep(retry_wait)
                retries += 1
            else:
                return f"Error: HTTP Error for {accession}: {str(err)}"
        except Exception as e:
            return f"Error: Entrez exception for {accession}: {str(e)}"
    return f"Error: Both methods failed for {accession}: Too Many Requests after retries."

def download_worker(accession, output_dir, email, api_key=None, resume=False):
    final_output_path = os.path.join(output_dir, f"{accession}.gbff")
    if resume and os.path.exists(final_output_path):
        return f"Skipped (already exists): {accession}"
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
                    final_output_path = os.path.join(output_dir, f"{accession}.gbff")
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

        # Fallback to Entrez with retry and rate limiting
        return fetch_entrez_genbank(accession, output_dir, email, api_key)

def main():
    args = parse_arguments()
    if not shutil.which("datasets"):
        print("Error: 'datasets' CLI tool not found in PATH. Please install NCBI Datasets CLI.")
        exit(1)
    os.makedirs(args.output, exist_ok=True)
    print(f"Reading accessions from {args.input}...")
    accessions = get_accessions_from_fasta(args.input)
    if not accessions:
        print("No valid accessions found. Exiting.")
        exit(0)
    api_key = args.api_key if args.api_key else None
    print(f"Found {len(accessions)} unique accessions. Starting parallel download with {args.threads} threads...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_acc = {
            executor.submit(download_worker, acc, args.output, args.email, api_key, args.resume): acc
            for acc in accessions
        }
        for future in as_completed(future_to_acc):
            result = future.result()
            print(result)
    print(f"\nProcessing complete. Files saved to: {args.output}")

if __name__ == "__main__":
    main()
