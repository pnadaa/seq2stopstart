#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import subprocess
import csv
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from typing import Dict, Any, List


def extract_sequence_with_blastdbcmd(accession: str, start: int, end: int, blast_db: str, temp_dir: str, verbose: bool = False) -> str:
    out_fasta = os.path.join(temp_dir, f"{accession}_{start}_{end}.fa")
    blast_cmd = [
        "blastdbcmd", "-db", blast_db,      # blast_db is a prefix e.g. /data/blastdb/mydb
        "-entry", accession,
        "-range", f"{start}-{end}",
        "-outfmt", "%f", "-out", out_fasta
    ]
    if verbose:
        print("Running: " + " ".join(blast_cmd))
    subprocess.run(blast_cmd, check=True)
    seq_record = next(SeqIO.parse(out_fasta, "fasta"))
    seq = str(seq_record.seq)
    os.remove(out_fasta)
    return seq


def align_query_to_ref(query_rna: str, ref_dna: str) -> Dict[str, Any] | None:
    """
    Local pairwise alignment of query against forward reference using
    Bio.Align.PairwiseAligner (replaces deprecated Bio.pairwise2).

    Scoring: match=2, mismatch=-1, gap open=-2, gap extend=-0.5
    Both query and its reverse complement are tried; the best-scoring
    alignment is returned. Alignment offsets are always into the forward
    reference sequence, so coordinate remapping downstream is correct
    regardless of IS orientation.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    query_dna = query_rna.replace('U', 'T')
    query_rc = str(Seq(query_dna).reverse_complement())

    best = None
    for strand, q in [('+', query_dna), ('-', query_rc)]:
        alignments = aligner.align(ref_dna, q)
        try:
            aln = next(iter(alignments))
        except StopIteration:
            continue
        # aln.aligned[0] holds (start, end) pairs for each aligned block
        # in the target (ref_dna); take outermost positions.
        ref_blocks = aln.aligned[0]
        s = {
            'strand': strand,
            'score': aln.score,
            'start': int(ref_blocks[0][0]),
            'end': int(ref_blocks[-1][1]),
            'aln': aln
        }
        if best is None or s['score'] > best['score']:
            best = s
    return best


def extract_genes_with_boundary(gb_record: Any, boundary_type: str) -> List[Dict[str, Any]]:
    genes = []
    for feat in gb_record.features:
        if feat.type in ('CDS', 'gene'):
            qual = feat.qualifiers
            gene = qual.get('locus_tag', ['?'])[0]
            strand = feat.location.strand
            loc = feat.location
            if boundary_type == "stop":
                boundary = loc.end - 1 if strand == 1 else loc.start
            else:
                boundary = loc.start if strand == 1 else loc.end - 1
            genes.append({
                'gene': gene,
                'start': int(loc.start),
                'end': int(loc.end),
                'strand': strand,
                'boundary': int(boundary)
            })
    return genes


def find_flanking_genes(query_boundary: int, genes: List[Dict[str, Any]]) -> tuple:
    min_upstream = None
    min_downstream = None
    min_dist_up = 1e12
    min_dist_down = 1e12
    for gene in genes:
        b = gene['boundary']
        if gene['strand'] == 1:
            dist_up = query_boundary - b
            if dist_up >= 0 and dist_up < min_dist_up:
                min_upstream, min_dist_up = gene, dist_up
            dist_down = b - query_boundary
            if dist_down >= 0 and dist_down < min_dist_down:
                min_downstream, min_dist_down = gene, dist_down
        elif gene['strand'] == -1:
            dist_up = b - query_boundary
            if dist_up >= 0 and dist_up < min_dist_up:
                min_upstream, min_dist_up = gene, dist_up
            dist_down = query_boundary - b
            if dist_down >= 0 and dist_down < min_dist_down:
                min_downstream, min_dist_down = gene, dist_down
    return (min_upstream, min_dist_up if min_upstream else None,
            min_downstream, min_dist_down if min_downstream else None)


def classify_alignment(align_start: int, align_end: int, genes: List[Dict[str, Any]]) -> str:
    """
    Returns one of:
      'inside'     – alignment fully contained within a single gene
      'partial'    – alignment partially overlaps one or more genes
      'intergenic' – no overlap with any gene

    Uses half-open interval arithmetic matching BioPython's SeqFeature
    coordinates. Completely independent of boundary_type.
    """
    overlapping = [
        g for g in genes
        if align_start < g['end'] and align_end > g['start']
    ]
    if not overlapping:
        return 'intergenic'
    for g in overlapping:
        if align_start >= g['start'] and align_end <= g['end']:
            return 'inside'
    return 'partial'


def _random_one_accession(
    accession: str,
    lengths: List[int],
    genbank_dir: str,
    gb_naming: str,
    boundary_type: str,
    n_random: int,
    seed: int,
    verbose: bool = False
) -> tuple[List[float], List[float]]:
    """
    Randomly place n_random sequences of sampled lengths within one genome.
    Each accession receives a deterministic but unique seed derived from the
    global seed so results are reproducible regardless of process order.
    Applies the same inside-gene and unannotated filters as the real pipeline.
    """
    rng = np.random.default_rng(seed)

    gb_path = (
        Path(genbank_dir) / f"{accession}.gbff"
        if gb_naming == "accession"
        else Path(genbank_dir) / gb_naming.format(accession=accession)
    )
    if not gb_path.exists():
        print(f"Warning: GenBank file not found for {accession}, skipping random control.")
        return [], []

    record = next(SeqIO.parse(gb_path, "genbank"))
    genome_length = len(record.seq)
    genes = extract_genes_with_boundary(record, boundary_type)

    if verbose:
        print(f"Generating {n_random} random placements for {accession} "
              f"(genome length: {genome_length:,} bp)")

    rand_ups, rand_downs = [], []
    placed = 0
    attempts = 0
    max_attempts = n_random * 20

    while placed < n_random and attempts < max_attempts:
        attempts += 1
        length = int(rng.choice(lengths))
        if length >= genome_length:
            continue
        rand_start = int(rng.integers(0, genome_length - length))
        rand_end = rand_start + length

        # Apply same geometric filter as real pipeline
        location = classify_alignment(rand_start, rand_end, genes)
        if location == 'inside':
            continue

        up_s, up_d_s, down_s, down_d_s = find_flanking_genes(rand_start, genes)
        up_e, up_d_e, down_e, down_d_e = find_flanking_genes(rand_end, genes)
        distances = [
            ("start", rand_start, up_s, up_d_s, down_s, down_d_s),
            ("end",   rand_end,   up_e, up_d_e, down_e, down_d_e)
        ]
        min_tuple = min(distances, key=lambda x: min(x[3] if x[2] else 1e12, x[5] if x[4] else 1e12))
        _, _, up_gene, up_dist, down_gene, down_dist = min_tuple

        # Apply same unannotated filter as real pipeline
        if up_dist is None and down_dist is None:
            continue

        if up_dist is not None:
            rand_ups.append(float(up_dist))
        if down_dist is not None:
            rand_downs.append(float(down_dist))
        placed += 1

    if placed < n_random:
        print(f"Warning: only placed {placed}/{n_random} random sequences for "
              f"{accession} after {max_attempts} attempts.")

    return rand_ups, rand_downs


def run_random_controls(
    results_filtered: List[Dict[str, Any]],
    genbank_dir: str,
    gb_naming: str,
    boundary_type: str,
    n_random: int,
    seed: int,
    nproc: int = 1,
    verbose: bool = False
) -> tuple[List[float], List[float]]:
    """
    Collects alignment lengths per accession from filtered real results, then
    dispatches per-accession random placement to _random_one_accession.
    Parallelised across accessions using the same --nproc as the real pipeline.
    Each accession gets a unique deterministic seed (seed + accession_index).
    """
    accession_lengths: Dict[str, List[int]] = defaultdict(list)
    for r in results_filtered:
        if not r.get('error') and r.get('align_start') is not None and r.get('align_end') is not None:
            length = r['align_end'] - r['align_start']
            if length > 0:
                accession_lengths[r['accession']].append(length)

    # Build task list with per-accession seeds: seed+i ensures reproducibility
    # regardless of which worker handles which accession.
    tasks = [
        (acc, lengths, genbank_dir, gb_naming, boundary_type, n_random, seed + i, verbose)
        for i, (acc, lengths) in enumerate(accession_lengths.items())
    ]

    if nproc > 1:
        with Pool(nproc) as pool:
            per_accession = pool.starmap(_random_one_accession, tasks)
    else:
        per_accession = [_random_one_accession(*t) for t in tasks]

    rand_ups   = [v for ups, _    in per_accession for v in ups]
    rand_downs = [v for _,  downs in per_accession for v in downs]
    return rand_ups, rand_downs


def process_one_region(region: Dict[str, Any], genbank_dir: str, blast_db: str, temp_dir: str, gb_naming: str, boundary_type: str, verbose: bool = False) -> Dict[str, Any]:
    try:
        accession = region['accession']
        start, end = region['start'], region['end']   # always start <= end after parse_fasta_regions
        is_reverse = region.get('is_reverse', False)

        # blastdbcmd always receives start <= end (normalised in parse_fasta_regions).
        # align_query_to_ref tries both query_dna and query_rc against the forward
        # reference, so reverse-strand entries align correctly without RC'ing ref_seq.
        # Alignment offsets are indices into the forward reference → remapping is correct.
        ref_seq = extract_sequence_with_blastdbcmd(accession, start, end, blast_db, temp_dir, verbose)
        best_aln = align_query_to_ref(region['sequence'], ref_seq)
        if not best_aln:
            raise RuntimeError("No good alignment found")

        coord_window_start = start          # normalised min coordinate
        rel_start = coord_window_start + best_aln['start']
        rel_end   = coord_window_start + best_aln['end']

        gb_path = (
            Path(genbank_dir) / f"{accession}.gbff"
            if gb_naming == "accession"
            else Path(genbank_dir) / gb_naming.format(accession=accession)
        )
        if not gb_path.exists():
            raise FileNotFoundError(f"Missing GenBank file: {gb_path}")
        record = next(SeqIO.parse(gb_path, "genbank"))
        genes = extract_genes_with_boundary(record, boundary_type)

        # Geometric classification: independent of boundary_type, so location
        # counts are consistent across --boundary_type stop and start runs.
        location = classify_alignment(rel_start, rel_end, genes)

        up_start, up_dist_start, down_start, down_dist_start = find_flanking_genes(rel_start, genes)
        up_end, up_dist_end, down_end, down_dist_end = find_flanking_genes(rel_end, genes)
        distances = [
            ("start", rel_start, up_start, up_dist_start, down_start, down_dist_start),
            ("end",   rel_end,   up_end,   up_dist_end,   down_end,   down_dist_end)
        ]
        min_tuple = min(distances, key=lambda x: min(x[3] if x[2] else 1e12, x[5] if x[4] else 1e12))
        which, rel_boundary, up_gene, up_dist, down_gene, down_dist = min_tuple

        up_gene_name       = up_gene['gene']      if up_gene   else None
        up_gene_boundary   = up_gene['boundary']  if up_gene   else None
        down_gene_name     = down_gene['gene']     if down_gene else None
        down_gene_boundary = down_gene['boundary'] if down_gene else None

        return {
            "accession": accession,
            "query": region.get('header'),
            "is_reverse": is_reverse,
            "location": location,
            "seq_start": start,
            "seq_end": end,
            "align_start": rel_start,
            "align_end": rel_end,
            "which_boundary_used": which,
            "boundary_used": rel_boundary,
            "up_gene": up_gene_name,
            "up_boundary": up_gene_boundary,
            "up_dist": up_dist,
            "down_gene": down_gene_name,
            "down_boundary": down_gene_boundary,
            "down_dist": down_dist,
            "score": best_aln['score'],
            "strand": best_aln['strand'],
            "error": None
        }
    except Exception as e:
        region_base = region if isinstance(region, dict) else {}
        return {
            "accession": region_base.get("accession"),
            "query": region_base.get("header"),
            "is_reverse": region_base.get("is_reverse"),
            "location": None,
            "seq_start": region_base.get("start"),
            "seq_end": region_base.get("end"),
            "align_start": None,
            "align_end": None,
            "which_boundary_used": None,
            "boundary_used": None,
            "up_gene": None,
            "up_boundary": None,
            "up_dist": None,
            "down_gene": None,
            "down_boundary": None,
            "down_dist": None,
            "score": None,
            "strand": None,
            "error": str(e)
        }


def _proc(region: Dict[str, Any], genbank_dir: str, blast_db: str, temp_dir: str, gb_naming: str, boundary_type: str, verbose: bool) -> Dict[str, Any]:
    return process_one_region(region, genbank_dir, blast_db, temp_dir, gb_naming, boundary_type, verbose)


def parse_fasta_regions(fasta_path: str) -> List[Dict[str, Any]]:
    results = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        header = rec.id
        try:
            parts = header.split('_')
            if len(parts) != 2:
                raise ValueError(f"Header does not have expected format (ACCESSION_START-END): {header}")
            coords = parts[1].split('-')
            if len(coords) != 2:
                raise ValueError(f"Coordinates not found or badly formatted in FASTA header: {header}")
            accession = parts[0]
            raw_start, raw_end = int(coords[0]), int(coords[1])

            # Detect reverse-strand entries (start > end) and normalise so that
            # blastdbcmd always receives a valid ascending range, and alignment
            # offsets remain indices into the forward sequence for correct remapping.
            is_reverse = raw_start > raw_end
            start = min(raw_start, raw_end)
            end   = max(raw_start, raw_end)

            results.append({
                "accession": accession,
                "start": start,
                "end": end,
                "is_reverse": is_reverse,
                "sequence": str(rec.seq),
                "header": header
            })
        except Exception as e:
            print(f"Error parsing FASTA header: {header}. Exception: {e}")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Region/genbank/BLAST pipeline with gene boundaries (start/stop) switchable by argument."
    )
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--genbank_dir", required=True)
    parser.add_argument(
        "--blast_db", required=True,
        help="Path to the BLAST+ database prefix (e.g. /data/blastdb/mydb), "
             "as passed to makeblastdb -out. Not a directory."
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--gb_naming", default="accession")
    parser.add_argument("--boundary_type", choices=["stop", "start"], default="stop",
                        help="Compare distance to either stop (default) or start codons of flanking genes")
    parser.add_argument(
        "--nproc", default=1, type=int,
        help="Number of parallel worker processes. On PBS, set this to match your "
             "ncpus allocation (e.g. #PBS -l ncpus=8 → --nproc 8). "
             "Applied to both real region processing and random control generation."
    )
    parser.add_argument("--csv_out", default="coordinates_with_genes.csv")
    parser.add_argument("--dist_out", default="distances.csv")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument(
        "--plot_name", default="distance_distribution.png",
        help="Filename for the distance distribution plot. "
             "The location category plot is saved with '_location' appended to the stem."
    )
    parser.add_argument(
        "--xlim", type=float, nargs=2, metavar=("XMIN", "XMAX"), default=None,
        help="Limit the x-axis of the distance plot, e.g. --xlim 0 500. "
             "Bar widths and ticks scale to this range."
    )
    parser.add_argument(
        "--n_random", type=int, default=1000,
        help="Number of random sequence placements per genome for the null distribution "
             "control. Lengths are sampled from real alignment lengths for that genome. "
             "Set to 0 to disable. Default: 1000."
    )
    parser.add_argument(
        "--random_seed", type=int, default=42,
        help="Random seed for reproducible control placement. Each accession gets a "
             "unique deterministic seed derived from this value. Default: 42."
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    print("\nArguments used:")
    for arg, val in vars(args).items():
        print(f"  {arg}: {val}")

    output_dir = Path(args.output_dir)
    temp_dir = output_dir / "blast_tmp"
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    if args.verbose:
        print(f"Temp dir: {temp_dir.resolve()}")

    regions = parse_fasta_regions(args.fasta)
    func = partial(
        _proc,
        genbank_dir=args.genbank_dir,
        blast_db=args.blast_db,
        temp_dir=str(temp_dir),
        gb_naming=args.gb_naming,
        boundary_type=args.boundary_type,
        verbose=args.verbose
    )
    if args.nproc > 1:
        with Pool(args.nproc) as pool:
            results = pool.map(func, regions)
    else:
        results = [func(region) for region in regions]

    # Capture counts from all valid (non-error) results before any filtering.
    # location classification is boundary_type-independent (geometric overlap).
    valid_before  = [r for r in results if not r['error']]
    n_inside      = sum(1 for r in valid_before if r.get('location') == 'inside')
    n_partial     = sum(1 for r in valid_before if r.get('location') == 'partial')
    n_unannotated = sum(1 for r in valid_before
                        if r.get('location') == 'intergenic'
                        and r.get('up_dist') is None
                        and r.get('down_dist') is None)
    n_outside     = len(valid_before) - n_inside - n_partial - n_unannotated

    # Write full CSV BEFORE filtering — all rows preserved with location values intact
    csv_out_path = output_dir / args.csv_out
    with open(csv_out_path, "w", newline="") as outcsv:
        fieldnames = [
            "accession", "query", "is_reverse", "location",
            "seq_start", "seq_end", "align_start", "align_end",
            "which_boundary_used", "boundary_used",
            "up_gene", "up_boundary", "up_dist",
            "down_gene", "down_boundary", "down_dist",
            "score", "strand", "error"
        ]
        writer = csv.DictWriter(outcsv, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    # Filter 1: alignment fully inside a gene body
    n_before = len(results)
    results_filtered = [r for r in results if r.get('location') != 'inside']
    n_removed = n_before - len(results_filtered)
    if n_removed:
        print(f"Filtered out {n_removed} result(s) where alignment is fully inside a gene body.")

    # Filter 2: unannotated regions (intergenic classification but no flanking genes, no error)
    n_before = len(results_filtered)
    results_filtered = [r for r in results_filtered
                        if r.get('error')
                        or not (r.get('up_dist') is None and r.get('down_dist') is None)]
    n_unannotated_removed = n_before - len(results_filtered)
    if n_unannotated_removed:
        print(f"Filtered out {n_unannotated_removed} result(s) with no flanking gene annotation.")

    # dist CSV uses filtered results (intergenic + partial only)
    dist_out_path = output_dir / args.dist_out
    with open(dist_out_path, "w", newline="") as dcsv:
        writer = csv.writer(dcsv)
        writer.writerow(["up_dist", "down_dist"])
        for r in results_filtered:
            if r['error']:
                continue
            writer.writerow([r['up_dist'], r['down_dist']])

    if args.plot:
        plot_stem   = Path(args.plot_name).stem
        plot_suffix = Path(args.plot_name).suffix or ".png"

        # --- Plot 1: distance distribution with random control overlay ---
        ups   = [float(r['up_dist'])   for r in results_filtered if r.get('up_dist')   is not None and not r.get('error')]
        downs = [float(r['down_dist']) for r in results_filtered if r.get('down_dist') is not None and not r.get('error')]

        # Generate random controls, parallelised across accessions via --nproc
        rand_ups, rand_downs = [], []
        if args.n_random > 0:
            print(f"Generating random controls ({args.n_random} placements per genome, "
                  f"nproc={args.nproc})...")
            rand_ups, rand_downs = run_random_controls(
                results_filtered,
                args.genbank_dir,
                args.gb_naming,
                args.boundary_type,
                args.n_random,
                args.random_seed,
                nproc=args.nproc,
                verbose=args.verbose
            )

        if not ups and not downs:
            print("Warning: no valid upstream/downstream distances to plot. Skipping distance plot.")
        else:
            fig1, ax1 = plt.subplots()

            # Use density when random controls are present so distributions with
            # different sample sizes are directly comparable on the same y-axis.
            use_density = bool(rand_ups or rand_downs)

            if args.xlim is not None:
                xmin, xmax = float(args.xlim[0]), float(args.xlim[1])
                n_bins = min(40, int(xmax - xmin))
                hist_kwargs = {'bins': n_bins, 'range': (xmin, xmax)}
            else:
                all_vals = ups + downs + rand_ups + rand_downs
                n_bins = min(40, int(max(all_vals) - min(all_vals)))
                hist_kwargs = {'bins': n_bins}

            if ups:
                ax1.hist(ups,   **hist_kwargs, alpha=0.7, label="Upstream (real)",
                         density=use_density)
            if downs:
                ax1.hist(downs, **hist_kwargs, alpha=0.7, label="Downstream (real)",
                         density=use_density)
            # Random controls as step outlines so they don't obscure real bars
            if rand_ups:
                ax1.hist(rand_ups,   **hist_kwargs, alpha=0.9, label="Upstream (random)",
                         density=use_density, histtype='step', linewidth=1.5, linestyle='--')
            if rand_downs:
                ax1.hist(rand_downs, **hist_kwargs, alpha=0.9, label="Downstream (random)",
                         density=use_density, histtype='step', linewidth=1.5, linestyle='--')

            if args.xlim is not None:
                ax1.set_xlim(xmin, xmax)
                ax1.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))

            ax1.set_xlabel(f"Distance to {args.boundary_type} codon (bp)")
            ax1.set_ylabel("Density" if use_density else "Count")
            ax1.set_title(f"{args.boundary_type.capitalize()} distance distribution "
                          f"({os.path.basename(args.fasta)})")
            ax1.legend()
            fig1.tight_layout()
            fig1.savefig(output_dir / args.plot_name)
            plt.close(fig1)
            if args.verbose:
                print(f"Distance plot saved: {(output_dir / args.plot_name).resolve()}")

        # --- Plot 2: location category counts (all valid results, pre-filter) ---
        fig2, ax2 = plt.subplots()
        labels = ['Inside gene', 'Partial overlap', 'Intergenic', 'Unannotated']
        counts = [n_inside, n_partial, n_outside, n_unannotated]
        bars = ax2.bar(labels, counts, alpha=0.7)
        y_offset = max(counts) * 0.01 if max(counts) > 0 else 0.5
        for bar, count in zip(bars, counts):
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + y_offset,
                str(count),
                ha='center', va='bottom'
            )
        ax2.set_ylabel("Count")
        ax2.set_title(f"Alignment location ({os.path.basename(args.fasta)})")
        fig2.tight_layout()
        loc_plot_name = f"{plot_stem}_location{plot_suffix}"
        fig2.savefig(output_dir / loc_plot_name)
        plt.close(fig2)
        if args.verbose:
            print(f"Location plot saved: {(output_dir / loc_plot_name).resolve()}")


if __name__ == "__main__":
    main()
