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
import tempfile
import csv
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from typing import Dict, Any, List


def extract_sequence_with_blastdbcmd(accession: str, start: int, end: int, blast_db: str, temp_dir: str, verbose: bool = False) -> str:
    # A unique temp file per call, not per coordinate window. A forward entry
    # (X_100-159) and its reverse counterpart (X_159-100) normalise to the same
    # window, so a name built from accession+coords collides between workers:
    # one would delete the file while another was still reading it, producing
    # sporadic FileNotFoundError / empty-record failures that differed from run
    # to run. mkstemp keeps each extraction private.
    fd, out_fasta = tempfile.mkstemp(prefix=f"{accession}_{start}_{end}_",
                                     suffix=".fa", dir=temp_dir)
    os.close(fd)
    try:
        blast_cmd = [
            "blastdbcmd", "-db", blast_db,      # blast_db is a prefix e.g. /data/blastdb/mydb
            "-entry", accession,
            "-range", f"{start}-{end}",
            "-outfmt", "%f", "-out", out_fasta
        ]
        if verbose:
            print("Running: " + " ".join(blast_cmd))
        subprocess.run(blast_cmd, check=True)
        seq_record = next(SeqIO.parse(out_fasta, "fasta"), None)
        if seq_record is None:
            # StopIteration would otherwise propagate with an empty message and
            # land in the CSV as a blank 'error', i.e. a silently dropped row.
            raise RuntimeError(
                f"blastdbcmd returned no sequence for {accession}:{start}-{end}"
            )
        return str(seq_record.seq)
    finally:
        if os.path.exists(out_fasta):
            os.remove(out_fasta)


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


def map_query_pos_to_ref(best_aln: Dict[str, Any], query_len: int, anchor_pos: int) -> int | None:
    """
    Map a 1-based position within the query (in the query's own, as-supplied
    orientation) onto a 0-based offset into the forward reference window.

    align_query_to_ref aligns either the query or its reverse complement against
    the forward reference; when the reverse complement won, position p of the
    original query sits at offset (query_len - 1 - (p - 1)) of the aligned
    sequence, so the anchor is remapped before the block lookup.

    Returns None when the anchor base is not covered by the local alignment
    (gapped, or outside the aligned span) — the caller falls back and flags it.
    """
    q0 = anchor_pos - 1
    if not (0 <= q0 < query_len):
        return None
    if best_aln['strand'] == '-':
        q0 = query_len - 1 - q0

    aln = best_aln['aln']
    for (t_start, _t_end), (q_start, q_end) in zip(aln.aligned[0], aln.aligned[1]):
        if q_start <= q0 < q_end:
            return int(t_start + (q0 - q_start))
    return None


def _gene_segments(parts: List[Any], strand: int | None) -> List[tuple]:
    """
    Half-open genomic intervals actually occupied by a feature.

    A feature that crosses the origin of a circular molecule, e.g.
    join(4044092..4044757,1..936), has BioPython .start == 0 and .end == the
    genome length, so treating it as one [start, end) interval makes that single
    gene "contain" the whole chromosome: every target in the genome was then
    classified 'inside' and dropped. Parts are listed in biological order, so the
    origin is crossed wherever consecutive parts step backwards (leftwards on the
    + strand, rightwards on the - strand); each run between such steps is one
    contiguous interval. Ordinary multi-part features (ribosomal frameshifts)
    never step backwards and stay a single interval.
    """
    runs: List[List[Any]] = [[parts[0]]]
    for prev, cur in zip(parts, parts[1:]):
        backwards = cur.start < prev.start if strand != -1 else cur.start > prev.start
        if backwards:
            runs.append([cur])
        else:
            runs[-1].append(cur)
    return [(int(min(p.start for p in run)), int(max(p.end for p in run))) for run in runs]


def extract_genes_with_boundary(gb_record: Any, boundary_type: str) -> List[Dict[str, Any]]:
    genes = []
    for feat in gb_record.features:
        if feat.type in ('CDS', 'gene') and feat.location is not None:
            qual = feat.qualifiers
            gene = qual.get('locus_tag', ['?'])[0]
            strand = feat.location.strand
            # Codons come from the first/last part in biological order, not from
            # the feature's min/max extent, which for an origin-spanning gene are
            # the two ends of the genome rather than its start and stop codons.
            parts = feat.location.parts
            first, last = parts[0], parts[-1]
            if boundary_type == "stop":
                boundary = last.end - 1 if strand == 1 else last.start
            else:
                boundary = first.start if strand == 1 else first.end - 1
            genes.append({
                'gene': gene,
                'start': int(feat.location.start),
                'end': int(feat.location.end),
                'segments': _gene_segments(parts, strand),
                'strand': strand,
                'boundary': int(boundary)
            })
    return genes


def find_flanking_genes(query_boundary: int, genes: List[Dict[str, Any]],
                        genome_length: int | None = None, circular: bool = False) -> tuple:
    """
    Nearest gene boundary upstream and downstream of query_boundary, each
    judged in that gene's own orientation. Every gene is a candidate regardless
    of where the query sits: a query inside a gene is credited to whichever
    boundary is nearest on each side, whether that belongs to the gene it is
    inside or to a neighbour.

    On a circular molecule distances are taken around the origin, so a query
    near either end of the sequence still reaches the next gene across it
    instead of reporting a spurious long distance or no gene at all.
    """
    wrap = circular and genome_length
    min_upstream = None
    min_downstream = None
    min_dist_up = 1e12
    min_dist_down = 1e12
    for gene in genes:
        b = gene['boundary']
        if gene['strand'] == 1:
            dist_up = query_boundary - b
            dist_down = b - query_boundary
        elif gene['strand'] == -1:
            dist_up = b - query_boundary
            dist_down = query_boundary - b
        else:
            continue
        if wrap:
            dist_up %= genome_length
            dist_down %= genome_length
        if dist_up >= 0 and dist_up < min_dist_up:
            min_upstream, min_dist_up = gene, dist_up
        if dist_down >= 0 and dist_down < min_dist_down:
            min_downstream, min_dist_down = gene, dist_down
    return (min_upstream, min_dist_up if min_upstream else None,
            min_downstream, min_dist_down if min_downstream else None)


def gene_containing(pos: int, genes: List[Dict[str, Any]]) -> str | None:
    """Locus tag(s) of every gene whose body covers 0-based position pos."""
    tags = sorted({g['gene'] for g in genes
                   if any(s <= pos < e for s, e in g['segments'])})
    return ";".join(tags) if tags else None


def classify_alignment(align_start: int, align_end: int, genes: List[Dict[str, Any]]) -> str:
    """
    Returns one of:
      'inside'     – alignment fully contained within a single gene
      'partial'    – alignment partially overlaps one or more genes
      'intergenic' – no overlap with any gene

    Uses half-open interval arithmetic matching BioPython's SeqFeature
    coordinates, against each gene's real segments (see _gene_segments).
    Completely independent of boundary_type, and purely descriptive: it does not
    decide whether a target is counted.
    """
    overlapping = [
        (s, e) for g in genes for s, e in g['segments']
        if align_start < e and align_end > s
    ]
    if not overlapping:
        return 'intergenic'
    for s, e in overlapping:
        if align_start >= s and align_end <= e:
            return 'inside'
    return 'partial'


def load_genbank(gb_path: Path, boundary_type: str) -> tuple[List[Dict[str, Any]], int, bool]:
    """Genes, sequence length and circular topology for one GenBank record."""
    record = next(SeqIO.parse(gb_path, "genbank"))
    circular = record.annotations.get('topology', '').lower() == 'circular'
    return extract_genes_with_boundary(record, boundary_type), len(record.seq), circular


def _random_one_accession(
    accession: str,
    lengths: List[int],
    genbank_dir: str,
    gb_naming: str,
    boundary_type: str,
    n_random: int,
    seed: int,
    verbose: bool = False,
    anchor_mode: str = "endpoints",
    anchor_pos: int = 31
) -> tuple[List[tuple], Dict[str, int]]:
    """
    Randomly place n_random sequences of sampled lengths within one genome.
    Each accession receives a deterministic but unique seed derived from the
    global seed so results are reproducible regardless of process order.
    Applies the same filter as the real pipeline: placements are only rejected
    when no flanking gene exists at all. Placements inside a gene are kept, as
    real targets inside a gene are.

    The measurement point matches the real pipeline exactly: with
    anchor_mode='center' distances are taken from a single anchor at the same
    offset into each random placement that the real anchor takes into the real
    target, rather than from whichever placement endpoint happens to be closer.
    Without this the null and the real distribution would not be comparable.

    Returns (rand_pairs, counts). rand_pairs holds one
    (accession, up_dist, down_dist) tuple per accepted placement — kept paired,
    and possibly containing None, so the written null CSV has the same row
    semantics as distances.csv. The accession label lets the null be collapsed
    per genome the same way the real data can be, which is needed to guard
    against pseudoreplication across targets from the same genome.
    counts tracks every valid placement attempt by location category for Plot 2b.
      counts['inside']     – accepted: fully inside a gene
      counts['partial']    – accepted: partially overlapped a gene
      counts['intergenic'] – accepted: fully intergenic and annotated
      counts['unannotated']– rejected: no flanking gene found in either direction
    """
    rng = np.random.default_rng(seed)

    gb_path = (
        Path(genbank_dir) / f"{accession}.gbff"
        if gb_naming == "accession"
        else Path(genbank_dir) / gb_naming.format(accession=accession)
    )
    if not gb_path.exists():
        print(f"Warning: GenBank file not found for {accession}, skipping random control.")
        return [], {'inside': 0, 'partial': 0, 'intergenic': 0, 'unannotated': 0}

    genes, genome_length, circular = load_genbank(gb_path, boundary_type)

    if verbose:
        print(f"Generating {n_random} random placements for {accession} "
              f"(genome length: {genome_length:,} bp)")

    rand_pairs: List[tuple] = []
    counts = {'inside': 0, 'partial': 0, 'intergenic': 0, 'unannotated': 0}
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

        location = classify_alignment(rand_start, rand_end, genes)

        if anchor_mode == "center":
            off = anchor_pos - 1 if 0 <= anchor_pos - 1 < length else length // 2
            up_gene, up_dist, down_gene, down_dist = find_flanking_genes(
                rand_start + off, genes, genome_length, circular)
        else:
            up_s, up_d_s, down_s, down_d_s = find_flanking_genes(rand_start, genes, genome_length, circular)
            up_e, up_d_e, down_e, down_d_e = find_flanking_genes(rand_end, genes, genome_length, circular)
            distances = [
                ("start", rand_start, up_s, up_d_s, down_s, down_d_s),
                ("end",   rand_end,   up_e, up_d_e, down_e, down_d_e)
            ]
            min_tuple = min(distances, key=lambda x: min(x[3] if x[2] else 1e12, x[5] if x[4] else 1e12))
            _, _, up_gene, up_dist, down_gene, down_dist = min_tuple

        if up_dist is None and down_dist is None:
            counts['unannotated'] += 1
            continue

        # Accepted placement — record location and distances
        counts[location] += 1   # 'inside', 'partial' or 'intergenic'
        rand_pairs.append((
            accession,
            float(up_dist) if up_dist is not None else None,
            float(down_dist) if down_dist is not None else None
        ))
        placed += 1

    if placed < n_random:
        print(f"Warning: only placed {placed}/{n_random} random sequences for "
              f"{accession} after {max_attempts} attempts.")

    return rand_pairs, counts


def run_random_controls(
    results_filtered: List[Dict[str, Any]],
    genbank_dir: str,
    gb_naming: str,
    boundary_type: str,
    n_random: int,
    seed: int,
    nproc: int = 1,
    verbose: bool = False,
    anchor_mode: str = "endpoints",
    anchor_pos: int = 31
) -> tuple[List[float], List[float], Dict[str, int], List[tuple]]:
    """
    Collects alignment lengths per accession from filtered real results, then
    dispatches per-accession random placement to _random_one_accession.
    Parallelised across accessions using the same --nproc as the real pipeline.
    Each accession gets a unique deterministic seed (seed + accession_index).
    Returns (rand_ups, rand_downs, aggregated_counts, rand_pairs).
    """
    accession_lengths: Dict[str, List[int]] = defaultdict(list)
    for r in results_filtered:
        if not r.get('error') and r.get('align_start') is not None and r.get('align_end') is not None:
            # align_start/align_end are 1-based inclusive on output, so a 60 bp
            # alignment spans end - start + 1 bases.
            length = r['align_end'] - r['align_start'] + 1
            if length > 0:
                accession_lengths[r['accession']].append(length)

    # Per-accession seeds: seed+i ensures reproducibility regardless of
    # which worker process handles which accession.
    tasks = [
        (acc, lengths, genbank_dir, gb_naming, boundary_type, n_random, seed + i,
         verbose, anchor_mode, anchor_pos)
        for i, (acc, lengths) in enumerate(accession_lengths.items())
    ]

    if nproc > 1:
        with Pool(nproc) as pool:
            per_accession = pool.starmap(_random_one_accession, tasks)
    else:
        per_accession = [_random_one_accession(*t) for t in tasks]

    rand_pairs = [p for pairs, _ in per_accession for p in pairs]
    rand_ups   = [u for _, u, _ in rand_pairs if u is not None]
    rand_downs = [d for _, _, d in rand_pairs if d is not None]

    agg_counts: Dict[str, int] = {'inside': 0, 'partial': 0, 'intergenic': 0, 'unannotated': 0}
    for _, counts in per_accession:
        for k in agg_counts:
            agg_counts[k] += counts[k]

    return rand_ups, rand_downs, agg_counts, rand_pairs


def process_one_region(region: Dict[str, Any], genbank_dir: str, blast_db: str, temp_dir: str, gb_naming: str, boundary_type: str, verbose: bool = False, anchor_mode: str = "endpoints", anchor_pos: int = 31) -> Dict[str, Any]:
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

        # blastdbcmd -range is 1-based inclusive, so ref_seq[0] is 0-based genome
        # index start-1. All arithmetic below is 0-based half-open, matching
        # BioPython's SeqFeature coordinates, so distances against gene boundaries
        # are exact. Coordinates are converted back to 1-based only on output.
        win0 = start - 1
        rel_start = win0 + best_aln['start']    # 0-based, inclusive
        rel_end   = win0 + best_aln['end']      # 0-based, half-open

        gb_path = (
            Path(genbank_dir) / f"{accession}.gbff"
            if gb_naming == "accession"
            else Path(genbank_dir) / gb_naming.format(accession=accession)
        )
        if not gb_path.exists():
            raise FileNotFoundError(f"Missing GenBank file: {gb_path}")
        genes, genome_length, circular = load_genbank(gb_path, boundary_type)

        # Geometric classification: independent of boundary_type, so location
        # counts are consistent across --boundary_type stop and start runs.
        location = classify_alignment(rel_start, rel_end, genes)

        if anchor_mode == "center":
            # Measure from a single, biologically meaningful point: the centre of
            # the trimmed target (the presumed insertion site), mapped through the
            # alignment so gaps and reverse-strand entries are handled exactly.
            off = map_query_pos_to_ref(best_aln, len(region['sequence']), anchor_pos)
            if off is not None:
                rel_boundary = win0 + off
                anchor_source = "aligned"
            else:
                rel_boundary = (rel_start + rel_end) // 2
                anchor_source = "midpoint_fallback"
            which = "center"
            up_gene, up_dist, down_gene, down_dist = find_flanking_genes(
                rel_boundary, genes, genome_length, circular)
        else:
            anchor_source = "endpoint_min"
            up_start, up_dist_start, down_start, down_dist_start = find_flanking_genes(
                rel_start, genes, genome_length, circular)
            up_end, up_dist_end, down_end, down_dist_end = find_flanking_genes(
                rel_end, genes, genome_length, circular)
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

        # Coordinates out as 1-based inclusive; distances stay as computed.
        return {
            "accession": accession,
            "query": region.get('header'),
            "is_reverse": is_reverse,
            "location": location,
            "seq_start": start,
            "seq_end": end,
            "align_start": rel_start + 1,
            "align_end": rel_end,
            "anchor": rel_boundary + 1,
            "anchor_source": anchor_source,
            "anchor_in_gene": gene_containing(rel_boundary, genes),
            "which_boundary_used": which,
            "boundary_used": rel_boundary + 1,
            "up_gene": up_gene_name,
            "up_boundary": up_gene_boundary + 1 if up_gene_boundary is not None else None,
            "up_dist": up_dist,
            "down_gene": down_gene_name,
            "down_boundary": down_gene_boundary + 1 if down_gene_boundary is not None else None,
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
            "anchor": None,
            "anchor_source": None,
            "anchor_in_gene": None,
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


def _proc(region: Dict[str, Any], genbank_dir: str, blast_db: str, temp_dir: str, gb_naming: str, boundary_type: str, verbose: bool, anchor_mode: str, anchor_pos: int) -> Dict[str, Any]:
    return process_one_region(region, genbank_dir, blast_db, temp_dir, gb_naming, boundary_type, verbose, anchor_mode, anchor_pos)


def parse_fasta_regions(fasta_path: str) -> List[Dict[str, Any]]:
    results = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        header = rec.id
        try:
            # Split on the LAST separator only: '_' occurs inside many accessions
            # (NC_055040, NM_001126745, ...), so a plain split() would reject them
            # and silently drop the record. Mirrors get_genbank.py's parsing.
            if ':' in header:
                accession, coord_part = header.rsplit(':', 1)
            elif '_' in header:
                accession, coord_part = header.rsplit('_', 1)
            else:
                raise ValueError(f"Header does not have expected format (ACCESSION_START-END): {header}")
            coords = coord_part.split('-')
            if len(coords) != 2:
                raise ValueError(f"Coordinates not found or badly formatted in FASTA header: {header}")
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
        "--anchor", choices=["endpoints", "center"], default="endpoints",
        help="Which point of the mapped target distances are measured from. "
             "'endpoints' (default, legacy behaviour): try both alignment ends and "
             "report whichever is closest to a gene boundary. 'center': measure from "
             "a single position inside the target (see --anchor_pos), i.e. the "
             "presumed insertion site. 'center' is the unbiased choice — 'endpoints' "
             "is a min-of-two statistic that shifts distances downward."
    )
    parser.add_argument(
        "--anchor_pos", type=int, default=31,
        help="1-based position within the trimmed target used as the measurement "
             "point when --anchor center, interpreted in the target's own "
             "(stranded) orientation and mapped through the alignment. Default 31, "
             "the first base past the midpoint of a 60 bp target."
    )
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
        help="Filename stem for distance distribution plots. "
             "Real data: {stem}.png. Random control: {stem}_random.png. "
             "Overlay (--overlay): {stem}_overlay.png. "
             "Location plots: {stem}_location.png and {stem}_location_random.png."
    )
    parser.add_argument(
        "--xlim", type=float, nargs=2, metavar=("XMIN", "XMAX"), default=None,
        help="Limit the x-axis of all distance plots, e.g. --xlim 0 500. "
             "Bar widths and ticks scale to this range. Applied independently "
             "to each separate plot, and to the overlay when --overlay is set."
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
    parser.add_argument(
        "--overlay", action="store_true",
        help="When set, produce a single overlay plot of real vs random distributions "
             "using density on the y-axis and step outlines for the random control, "
             "instead of two separate count-based plots. "
             "Saves as {stem}_overlay.png instead of {stem}.png + {stem}_random.png."
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
        verbose=args.verbose,
        anchor_mode=args.anchor,
        anchor_pos=args.anchor_pos
    )
    if args.nproc > 1:
        with Pool(args.nproc) as pool:
            results = pool.map(func, regions)
    else:
        results = [func(region) for region in regions]

    # Capture counts from all valid (non-error) results before any filtering.
    # location classification is boundary_type-independent (geometric overlap).
    # A row with no flanking gene on either side has nothing to measure against
    # (in practice a GenBank record with no gene/CDS features), whatever its
    # location label; every other row is counted, inside a gene or not.
    valid_before  = [r for r in results if not r['error']]
    annotated     = [r for r in valid_before
                     if not (r.get('up_dist') is None and r.get('down_dist') is None)]
    n_unannotated = len(valid_before) - len(annotated)
    n_inside      = sum(1 for r in annotated if r.get('location') == 'inside')
    n_partial     = sum(1 for r in annotated if r.get('location') == 'partial')
    n_outside     = sum(1 for r in annotated if r.get('location') == 'intergenic')

    # Write full CSV BEFORE filtering — all rows preserved with location values intact
    csv_out_path = output_dir / args.csv_out
    with open(csv_out_path, "w", newline="") as outcsv:
        fieldnames = [
            "accession", "query", "is_reverse", "location",
            "seq_start", "seq_end", "align_start", "align_end",
            "anchor", "anchor_source", "anchor_in_gene",
            "which_boundary_used", "boundary_used",
            "up_gene", "up_boundary", "up_dist",
            "down_gene", "down_boundary", "down_dist",
            "score", "strand", "error"
        ]
        writer = csv.DictWriter(outcsv, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    # Targets inside a gene are NOT filtered. Dropping them removed exactly the
    # sites furthest from any codon and pulled the distance distribution towards
    # zero. Their distances run to the nearest boundary on each side, which may
    # belong to the gene they sit in or to a neighbouring gene.
    # Only filter: no flanking gene annotation at all (no error, both distances None).
    n_before = len(results)
    results_filtered = [r for r in results
                        if r.get('error')
                        or not (r.get('up_dist') is None and r.get('down_dist') is None)]
    n_unannotated_removed = n_before - len(results_filtered)
    if n_unannotated_removed:
        print(f"Filtered out {n_unannotated_removed} result(s) with no flanking gene annotation.")
    print(f"Distances reported for {n_inside + n_partial + n_outside} target(s): "
          f"{n_inside} inside a gene, {n_partial} partially overlapping, "
          f"{n_outside} intergenic.")

    # dist CSV uses filtered results (inside + partial + intergenic)
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

        ups   = [float(r['up_dist'])   for r in results_filtered if r.get('up_dist')   is not None and not r.get('error')]
        downs = [float(r['down_dist']) for r in results_filtered if r.get('down_dist') is not None and not r.get('error')]

        # Generate random controls, parallelised across accessions via --nproc
        rand_ups: List[float] = []
        rand_downs: List[float] = []
        rand_pairs: List[tuple] = []
        rand_counts: Dict[str, int] = {'inside': 0, 'partial': 0, 'intergenic': 0, 'unannotated': 0}
        if args.n_random > 0:
            print(f"Generating random controls ({args.n_random} placements per genome, "
                  f"nproc={args.nproc})...")
            rand_ups, rand_downs, rand_counts, rand_pairs = run_random_controls(
                results_filtered,
                args.genbank_dir,
                args.gb_naming,
                args.boundary_type,
                args.n_random,
                args.random_seed,
                nproc=args.nproc,
                verbose=args.verbose,
                anchor_mode=args.anchor,
                anchor_pos=args.anchor_pos
            )

            # Persist the null distances. Previously they were plotted and then
            # discarded, so any later test against the null needed a full re-run.
            # One row per accepted placement, same columns as distances.csv.
            rand_out_path = output_dir / f"{Path(args.dist_out).stem}_random{Path(args.dist_out).suffix or '.csv'}"
            with open(rand_out_path, "w", newline="") as rcsv:
                writer = csv.writer(rcsv)
                writer.writerow(["accession", "up_dist", "down_dist"])
                for acc, up_d, down_d in rand_pairs:
                    writer.writerow([
                        acc,
                        "" if up_d is None else up_d,
                        "" if down_d is None else down_d
                    ])
            print(f"Random control distances written: {rand_out_path.resolve()} "
                  f"({len(rand_pairs)} placements)")

        # --- Shared helpers (defined here so they close over args.xlim) ---

        def _make_hist_kwargs(vals: List[float]) -> dict:
            """
            Compute bins and range independently from the supplied data.
            When --xlim is set, all plots use the same visible window but
            each still derives its bin count from that window independently.
            When --xlim is not set, bins are computed from each dataset's
            own range, keeping the real and random plots fully independent.
            """
            if args.xlim is not None:
                xmin, xmax = float(args.xlim[0]), float(args.xlim[1])
                n_bins = min(100, int(xmax - xmin))
                return {'bins': n_bins, 'range': (xmin, xmax)}
            if not vals:
                return {'bins': 100}
            span = int(max(vals) - min(vals))
            return {'bins': min(100, span) if span > 0 else 1}

        def _apply_xlim_ticks(ax) -> None:
            """Apply xlim and scaled ticks only when --xlim is explicitly set."""
            if args.xlim is not None:
                ax.set_xlim(float(args.xlim[0]), float(args.xlim[1]))
                ax.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 5, 10]))

        def _annotate_bars(ax, bars, counts):
            y_offset = max(counts) * 0.01 if max(counts) > 0 else 0.5
            for bar, count in zip(bars, counts):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + y_offset,
                    str(count),
                    ha='center', va='bottom'
                )

        # ---------------------------------------------------------------
        # Distance plots — behaviour depends on --overlay
        # ---------------------------------------------------------------

        if args.overlay:
            # ---- Plot 1 (overlay): real + random on one figure, density y-axis ----
            # Bins computed from combined data so both distributions share the same grid.
            if not (ups or downs) and not (rand_ups or rand_downs):
                print("Warning: no distance data available. Skipping overlay plot.")
            else:
                combined = ups + downs + rand_ups + rand_downs
                hk_ov = _make_hist_kwargs(combined)

                fig_ov, ax_ov = plt.subplots()
                if ups:
                    ax_ov.hist(ups,   **hk_ov, alpha=0.7, density=True, label="Upstream (real)")
                if downs:
                    ax_ov.hist(downs, **hk_ov, alpha=0.7, density=True, label="Downstream (real)")
                if rand_ups:
                    ax_ov.hist(rand_ups,   **hk_ov, density=True, histtype='step',
                               linewidth=1.5, linestyle='--', label="Upstream (random)")
                if rand_downs:
                    ax_ov.hist(rand_downs, **hk_ov, density=True, histtype='step',
                               linewidth=1.5, linestyle='--', label="Downstream (random)")
                _apply_xlim_ticks(ax_ov)
                ax_ov.set_xlabel(f"Distance to {args.boundary_type} codon (bp)")
                ax_ov.set_ylabel("Density")
                ax_ov.set_title(f"{args.boundary_type.capitalize()} distance: real vs random "
                                f"({os.path.basename(args.fasta)})")
                ax_ov.legend()
                fig_ov.tight_layout()
                overlay_name = f"{plot_stem}_overlay{plot_suffix}"
                fig_ov.savefig(output_dir / overlay_name)
                plt.close(fig_ov)
                if args.verbose:
                    print(f"Overlay plot saved: {(output_dir / overlay_name).resolve()}")

        else:
            # ---- Plot 1a: real distance distribution (counts) ----
            # Bins derived solely from real data; completely independent of random data.
            if not ups and not downs:
                print("Warning: no valid upstream/downstream distances to plot. Skipping real distance plot.")
            else:
                hk_real = _make_hist_kwargs(ups + downs)
                fig1, ax1 = plt.subplots()
                if ups:
                    ax1.hist(ups,   **hk_real, alpha=0.7, label="Upstream")
                if downs:
                    ax1.hist(downs, **hk_real, alpha=0.7, label="Downstream")
                _apply_xlim_ticks(ax1)
                ax1.set_xlabel(f"Distance to {args.boundary_type} codon (bp)")
                ax1.set_ylabel("Count")
                ax1.set_title(f"{args.boundary_type.capitalize()} distance distribution "
                              f"({os.path.basename(args.fasta)})")
                ax1.legend()
                fig1.tight_layout()
                fig1.savefig(output_dir / args.plot_name)
                plt.close(fig1)
                if args.verbose:
                    print(f"Real distance plot saved: {(output_dir / args.plot_name).resolve()}")

            # ---- Plot 1b: random distance distribution (counts) ----
            # Bins derived solely from random data; completely independent of real data.
            if rand_ups or rand_downs:
                hk_rand = _make_hist_kwargs(rand_ups + rand_downs)
                fig_r, ax_r = plt.subplots()
                if rand_ups:
                    ax_r.hist(rand_ups,   **hk_rand, alpha=0.7, label="Upstream (random)")
                if rand_downs:
                    ax_r.hist(rand_downs, **hk_rand, alpha=0.7, label="Downstream (random)")
                _apply_xlim_ticks(ax_r)
                ax_r.set_xlabel(f"Distance to {args.boundary_type} codon (bp)")
                ax_r.set_ylabel("Count")
                ax_r.set_title(f"{args.boundary_type.capitalize()} random control distance distribution "
                               f"({os.path.basename(args.fasta)})")
                ax_r.legend()
                fig_r.tight_layout()
                rand_plot_name = f"{plot_stem}_random{plot_suffix}"
                fig_r.savefig(output_dir / rand_plot_name)
                plt.close(fig_r)
                if args.verbose:
                    print(f"Random distance plot saved: {(output_dir / rand_plot_name).resolve()}")

        # ---------------------------------------------------------------
        # Location plots — always separate, unaffected by --overlay
        # ---------------------------------------------------------------

        # ---- Plot 2a: real location category counts ----
        fig2, ax2 = plt.subplots()
        loc_labels = ['Inside gene', 'Partial overlap', 'Intergenic', 'Unannotated']
        real_counts = [n_inside, n_partial, n_outside, n_unannotated]
        bars2 = ax2.bar(loc_labels, real_counts, alpha=0.7)
        _annotate_bars(ax2, bars2, real_counts)
        ax2.set_ylabel("Count")
        ax2.set_title(f"Alignment location ({os.path.basename(args.fasta)})")
        fig2.tight_layout()
        loc_plot_name = f"{plot_stem}_location{plot_suffix}"
        fig2.savefig(output_dir / loc_plot_name)
        plt.close(fig2)
        if args.verbose:
            print(f"Real location plot saved: {(output_dir / loc_plot_name).resolve()}")

        # ---- Plot 2b: random location category counts ----
        # Shows the proportion of random placements falling in each category,
        # including placements rejected as unannotated, for comparison with the
        # real data.
        if any(rand_counts.values()):
            fig2r, ax2r = plt.subplots()
            rand_loc_counts = [
                rand_counts['inside'],
                rand_counts['partial'],
                rand_counts['intergenic'],
                rand_counts['unannotated']
            ]
            bars2r = ax2r.bar(loc_labels, rand_loc_counts, alpha=0.7)
            _annotate_bars(ax2r, bars2r, rand_loc_counts)
            ax2r.set_ylabel("Count")
            ax2r.set_title(f"Random control alignment location ({os.path.basename(args.fasta)})")
            fig2r.tight_layout()
            loc_rand_name = f"{plot_stem}_location_random{plot_suffix}"
            fig2r.savefig(output_dir / loc_rand_name)
            plt.close(fig2r)
            if args.verbose:
                print(f"Random location plot saved: {(output_dir / loc_rand_name).resolve()}")


if __name__ == "__main__":
    main()
