#!/usr/bin/env python3
"""
Aggregate the per-clade seq2startstop runs in one results directory.

Reads each <clade>_<boundary>/coordinates_with_genes.csv (real) and
distances_random.csv (matched null), and produces:

  clade_summary.csv   one row per clade x boundary_type x direction
  ks_tests.csv        significance tests; the "statistic" column is the KS D for
                      the KS families and the Wilcoxon W for the paired family
  prism_cdf_clades.csv        one column per clade x boundary x direction, every
                              analysed distance, for cumulative plots in Prism
  ecdf_families_x<N>.png      IS4 vs IS1182 ECDFs
  ecdf_main_clades_x<N>.png   ECDFs for the families and clades A/B/C_all
  ecdf_by_location_x<N>.png   IS1182_all / IS4_all split by where the target sits
  ecdf_cladeC_subclades_{start,stop}_x<N>.png    small multiples per subclade
  CAVEATS.md          interpretation caveats that travel with the numbers

Every target with a flanking gene annotation is analysed, including targets
that sit inside a gene; their distances run to the nearest boundary on each
side, whichever gene it belongs to. Distances are measured from the centre of
the 60 bp trimmed target (--anchor center --anchor_pos 31); the null uses the
same anchor.

Run under the biopytools environment (pandas + scipy):
  micromamba run -n biopytools python summarise_clades.py --results_dir results/IS1182_clades_v2
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# Categorical slots 1-5 of the validated default palette, in documented order.
# Validated on the adjacent pairlist (the one that applies to line charts):
# worst adjacent CVD dE 9.1, worst adjacent normal-vision dE 19.6 on the light
# surface. Slots 3-5 sit below 3:1 contrast, so the relief rule applies and the
# figures carry direct end-labels plus clade_summary.csv as the table view.
# Anything with more groups than this is drawn as small multiples, not more hues.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
NULL_INK = "#8a8a86"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID_INK = "#dcdcd8"

FAMILIES = ["IS4_all", "IS1182_all"]
MAIN_CLADES = ["cladeA", "cladeB", "cladeC_all"]
# Fixed colour per entity, so a series keeps its hue across every figure.
MAIN_SERIES = FAMILIES + MAIN_CLADES
SERIES_COLOUR = dict(zip(MAIN_SERIES, SERIES))

# Pairs that are genuinely disjoint and so may be compared directly.
# IS4_all vs IS1182_all are different IS families; cladeA/B/C_all partition
# IS1182_all exactly, so a clade may never be compared against IS1182_all.
DISJOINT_PAIRS = [
    ("IS4_all", "IS1182_all"),
    ("cladeA", "cladeB"), ("cladeA", "cladeC_all"), ("cladeB", "cladeC_all"),
]
SUBCLADES = [
    "cladeC_fba3", "cladeC_pcc4",
    "cladeC_unknown1", "cladeC_unknown2", "cladeC_unknown3",
    "cladeC_unknown4", "cladeC_unknown5",
]
BOUNDARIES = ["start", "stop"]
DIRECTIONS = ["up", "down"]
DIR_LABEL = {"up": "Upstream", "down": "Downstream"}

FASTA_STEM = "trimmed_remapped_60bp_trimmed"
FASTA_FOR = {
    "cladeA": f"{FASTA_STEM}_cladeA_targets.fasta",
    "cladeB": f"{FASTA_STEM}_cladeB_targets.fasta",
    "cladeC_all": f"{FASTA_STEM}_cladeC_all_targets.fasta",
    "IS1182_all": f"{FASTA_STEM}_IS1182_stranded_unique_targets_dedup.fasta",
    "IS4_all": f"{FASTA_STEM}_IS4_stranded_unique_targets_dedup.fasta",
    **{s: f"{FASTA_STEM}_{s}_targets.fasta" for s in SUBCLADES},
}


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def count_fasta_records(path: Path) -> int:
    if not path.exists():
        return -1
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def load_run(run_dir: Path) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Return (real, null) frames for one <clade>_<boundary> directory."""
    real_path = run_dir / "coordinates_with_genes.csv"
    null_path = run_dir / "distances_random.csv"
    real = pd.read_csv(real_path) if real_path.exists() else None
    null = pd.read_csv(null_path) if null_path.exists() else None
    return real, null


def analysed_subset(real: pd.DataFrame) -> pd.DataFrame:
    """
    Rows the pipeline actually reports distances for: no error and at least one
    flanking gene annotated. Targets inside a gene are kept. Mirrors the filter
    applied in seq2startstop.main before distances.csv is written.
    """
    ok = real[real["error"].isna()]
    return ok[~(ok["up_dist"].isna() & ok["down_dist"].isna())]


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------

def describe(values: pd.Series) -> dict:
    v = pd.to_numeric(values, errors="coerce").dropna()
    if len(v) == 0:
        return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan,
                    frac_le50=np.nan, frac_le100=np.nan, frac_le200=np.nan)
    return dict(
        n=int(len(v)),
        median=float(v.median()),
        q1=float(v.quantile(0.25)),
        q3=float(v.quantile(0.75)),
        frac_le50=float((v <= 50).mean()),
        frac_le100=float((v <= 100).mean()),
        frac_le200=float((v <= 200).mean()),
    )


def ks(a: pd.Series, b: pd.Series) -> tuple[float, float, int, int]:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan, len(a), len(b)
    res = stats.ks_2samp(a, b)
    return float(res.statistic), float(res.pvalue), len(a), len(b)


def paired_wilcoxon(real: pd.Series, null: pd.Series) -> tuple[float, float, int, float]:
    """
    Wilcoxon signed-rank on per-genome medians, matched by accession.

    The null is generated inside each genome, so real and null values are paired
    by genome rather than independent. A paired test therefore both respects the
    design and removes between-genome variation in gene density, which is the
    dominant nuisance term. Returns (statistic, p, n_pairs, median_difference)
    where the difference is real - null.
    """
    joined = pd.concat([real.rename("real"), null.rename("null")], axis=1,
                       join="inner").dropna()
    if len(joined) < 6:                       # Wilcoxon is meaningless below this
        return np.nan, np.nan, len(joined), np.nan
    diff = joined["real"] - joined["null"]
    if (diff == 0).all():
        return np.nan, np.nan, len(joined), 0.0
    res = stats.wilcoxon(joined["real"], joined["null"], zero_method="wilcox")
    return float(res.statistic), float(res.pvalue), len(joined), float(diff.median())


def paired_sign_test(real: pd.Series, null: pd.Series) -> tuple[float, float, int]:
    """
    Sign test on per-genome medians: in what fraction of genomes is the observed
    distance smaller than that genome's own null?

    Reported alongside the Wilcoxon because these distances are bounded below by
    zero and unbounded above, so a few genomes sitting very far from any codon
    dominate the magnitude-weighted signed ranks and can cancel out a consistent
    but modest shift in the majority. Where the two disagree, the sign test is
    describing how many genomes show the effect and the Wilcoxon how large the
    signed differences are; both are worth reporting.

    Returns (fraction_below_null, two_sided_p, n_pairs_excluding_ties).
    """
    joined = pd.concat([real.rename("real"), null.rename("null")], axis=1,
                       join="inner").dropna()
    diff = joined["real"] - joined["null"]
    nonzero = diff[diff != 0]
    n = len(nonzero)
    if n < 6:
        return np.nan, np.nan, n
    below = int((nonzero < 0).sum())
    p = stats.binomtest(below, n, 0.5, alternative="two-sided").pvalue
    return below / n, float(p), n


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni adjusted p-values; NaNs pass through untouched."""
    idx = [i for i, p in enumerate(pvals) if p == p]      # non-NaN
    out = [np.nan] * len(pvals)
    if not idx:
        return out
    order = sorted(idx, key=lambda i: pvals[i])
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        adj = min(1.0, (m - rank) * pvals[i])
        running = max(running, adj)                        # enforce monotonicity
        out[i] = running
    return out


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def ecdf(values) -> tuple[np.ndarray, np.ndarray]:
    v = np.sort(np.asarray(pd.to_numeric(pd.Series(values), errors="coerce").dropna(), dtype=float))
    if v.size == 0:
        return np.array([]), np.array([])
    return v, np.arange(1, v.size + 1) / v.size


def style_axis(ax, xmax):
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1.0)
    ax.grid(True, color=GRID_INK, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID_INK)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=8)


def place_end_labels(ax, entries, min_gap=0.058, lo=0.03, hi=0.97):
    """
    Draw right-edge series labels, nudged apart so converging curves do not
    stack their labels on top of each other. `entries` is [(y, text, colour)]
    in axes-fraction y.

    Labels are spread greedily upward from the lowest, then the whole stack is
    slid back down if it overflowed the top; if even that will not fit, the gap
    is shrunk to whatever the axis can hold. Without the slide-back step a
    crowded panel just re-collides at y=1.
    """
    if not entries:
        return
    entries = sorted(entries, key=lambda e: e[0])
    gap = min(min_gap, (hi - lo) / max(1, len(entries) - 1))

    placed = []
    for y, _, _ in entries:
        y = min(max(y, lo), hi)
        if placed and y - placed[-1] < gap:
            y = placed[-1] + gap
        placed.append(y)

    overflow = placed[-1] - hi
    if overflow > 0:                      # slide the stack down, keeping spacing
        shift = min(overflow, placed[0] - lo)
        placed = [y - shift for y in placed]
        placed = [min(max(y, lo), hi) for y in placed]

    for y, (_, text, colour) in zip(placed, entries):
        ax.annotate(text, xy=(1.0, y), xycoords="axes fraction",
                    xytext=(-4, 0), textcoords="offset points",
                    color=colour, fontsize=8, ha="right", va="center",
                    fontweight="bold")


def plot_series_grid(data, series, xmax, out_path, title, subtitle=None, ncol_legend=3):
    """
    2x2: rows = boundary type, cols = direction. One coloured line per entry in
    `series`, with its matched null dashed in the same hue. Used for both the
    family-level figure and the combined families+clades figure.
    """
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), sharex=True, sharey=True)
    fig.patch.set_facecolor("#fcfcfb")

    for r, boundary in enumerate(BOUNDARIES):
        for c, direction in enumerate(DIRECTIONS):
            ax = axes[r][c]
            ax.set_facecolor("#fcfcfb")
            end_labels = []
            for name in series:
                entry = data.get((name, boundary))
                if entry is None:
                    continue
                colour = SERIES_COLOUR[name]
                x, y = ecdf(entry["real"][f"{direction}_dist"])
                if x.size:
                    ax.plot(x, y, color=colour, linewidth=2, label=f"{name} (n={x.size})")
                    # Relief rule: direct end-labels, so identity never rests on
                    # colour alone (slots 3-5 sit under 3:1 on the light surface).
                    inside = x[x <= xmax]
                    if inside.size:
                        end_labels.append((float(y[inside.size - 1]), name, colour))
                if entry["null"] is not None:
                    xr, yr = ecdf(entry["null"][f"{direction}_dist"])
                    if xr.size:
                        ax.plot(xr, yr, color=colour, linewidth=1.1, linestyle="--",
                                alpha=0.5, label=f"{name} null")
            style_axis(ax, xmax)
            place_end_labels(ax, end_labels)
            if r == 0:
                ax.set_title(f"{DIR_LABEL[direction]}", color=TEXT_PRIMARY, fontsize=11)
            if c == 0:
                ax.set_ylabel(f"{boundary} codon\ncumulative fraction",
                              color=TEXT_SECONDARY, fontsize=9)
            if r == 1:
                ax.set_xlabel("Distance from target centre (bp)",
                              color=TEXT_SECONDARY, fontsize=9)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=ncol_legend, frameon=False,
               fontsize=8, labelcolor=TEXT_SECONDARY, bbox_to_anchor=(0.5, -0.005))
    full_title = title + "\nsolid = observed, dashed = matched random-placement null"
    if subtitle:
        full_title += f"\n{subtitle}"
    fig.suptitle(full_title, color=TEXT_PRIMARY, fontsize=12)
    fig.tight_layout(rect=(0, 0.08, 1, 0.93))
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out_path}")


LOCATIONS = [("all", "all targets"), ("inside", "inside a gene"),
             ("partial", "partial overlap"), ("intergenic", "intergenic")]


def plot_by_location(data, clade, xmax, out_path):
    """
    2x2 ECDFs for one set, one line per target location class plus all targets
    together, with the set's matched null. Shows where the targets that used to
    be dropped as 'inside a gene' fall relative to the rest.
    """
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), sharex=True, sharey=True)
    fig.patch.set_facecolor("#fcfcfb")

    for r, boundary in enumerate(BOUNDARIES):
        entry = data.get((clade, boundary))
        for c, direction in enumerate(DIRECTIONS):
            ax = axes[r][c]
            ax.set_facecolor("#fcfcfb")
            end_labels = []
            if entry is not None:
                col = f"{direction}_dist"
                real = entry["real"]
                for i, (key, label) in enumerate(LOCATIONS):
                    subset = real if key == "all" else real[real["location"] == key]
                    x, y = ecdf(subset[col])
                    if not x.size:
                        continue
                    ax.plot(x, y, color=SERIES[i], linewidth=2.4 if key == "all" else 1.6,
                            label=f"{label} (n={x.size})")
                    inside = x[x <= xmax]
                    if inside.size:
                        end_labels.append((float(y[inside.size - 1]), label, SERIES[i]))
                if entry["null"] is not None:
                    xr, yr = ecdf(entry["null"][col])
                    if xr.size:
                        ax.plot(xr, yr, color=NULL_INK, linewidth=1.1, linestyle="--",
                                label="random-placement null")
            style_axis(ax, xmax)
            place_end_labels(ax, end_labels)
            if r == 0:
                ax.set_title(DIR_LABEL[direction], color=TEXT_PRIMARY, fontsize=11)
            if c == 0:
                ax.set_ylabel(f"{boundary} codon\ncumulative fraction",
                              color=TEXT_SECONDARY, fontsize=9)
            if r == 1:
                ax.set_xlabel("Distance from target centre (bp)",
                              color=TEXT_SECONDARY, fontsize=9)

    # Legend from the panel with the most entries, so n reflects that panel.
    best = max((a for row in axes for a in row),
               key=lambda a: len(a.get_legend_handles_labels()[0]))
    handles, _ = best.get_legend_handles_labels()
    labels = [h.get_label().split(" (n=")[0] for h in handles]
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False,
               fontsize=8, labelcolor=TEXT_SECONDARY, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle(
        f"{clade}: distance from target centre to nearest flanking codon, by target location\n"
        "location = where the 60 bp target sits; every class is measured to the nearest "
        "boundary on each side, whichever gene owns it",
        color=TEXT_PRIMARY, fontsize=12,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out_path}")


def write_prism_cdf(data, out_path):
    """
    One column per <clade>_<Boundary>_<Direction>, holding every analysed
    distance for that run (blank-padded to the longest column), ready to paste
    into Prism for cumulative-distribution plots. Same rows as the ECDFs here.
    """
    order = FAMILIES + MAIN_CLADES + SUBCLADES
    columns = {}
    for clade in order + sorted({c for c, _ in data} - set(order)):
        for boundary in BOUNDARIES:
            entry = data.get((clade, boundary))
            if entry is None:
                continue
            for direction in DIRECTIONS:
                vals = pd.to_numeric(entry["real"][f"{direction}_dist"], errors="coerce").dropna()
                name = f"{clade}_{boundary.capitalize()}_{direction.capitalize()}"
                # Nullable ints: padding stays blank instead of turning every value into "12.0".
                columns[name] = vals.round().astype("Int64").reset_index(drop=True)
    pd.DataFrame(columns).to_csv(out_path, index=False)
    print(f"wrote {out_path}  ({len(columns)} columns)")


def plot_subclades(data, boundary, xmax, out_path):
    """Small multiples — one panel per cladeC subclade, plus cladeC_all."""
    panels = ["cladeC_all"] + SUBCLADES
    ncol = 4
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.2 * nrow),
                             sharex=True, sharey=True)
    fig.patch.set_facecolor("#fcfcfb")
    axes = np.atleast_1d(axes).ravel()

    for ax, clade in zip(axes, panels):
        ax.set_facecolor("#fcfcfb")
        entry = data.get((clade, boundary))
        if entry is None:
            ax.set_title(f"{clade}\n(no output)", color=TEXT_SECONDARY, fontsize=9)
            style_axis(ax, xmax)
            continue
        n_shown = 0
        for i, direction in enumerate(DIRECTIONS):
            x, y = ecdf(entry["real"][f"{direction}_dist"])
            if x.size:
                n_shown = max(n_shown, x.size)
                ax.plot(x, y, color=SERIES[i], linewidth=2, label=DIR_LABEL[direction])
            if entry["null"] is not None:
                xr, yr = ecdf(entry["null"][f"{direction}_dist"])
                if xr.size:
                    ax.plot(xr, yr, color=NULL_INK, linewidth=1.1,
                            linestyle="--" if direction == "up" else ":",
                            label=f"{DIR_LABEL[direction]} null")
        style_axis(ax, xmax)
        ax.set_title(f"{clade}  (n={n_shown})", color=TEXT_PRIMARY, fontsize=10)

    for ax in axes[len(panels):]:
        ax.set_visible(False)

    # x-label on the bottom row of visible panels only; a figure-level supxlabel
    # would land on top of the legend.
    for ax in axes[max(0, len(panels) - ncol):len(panels)]:
        ax.set_xlabel("Distance from target centre (bp)",
                      color=TEXT_SECONDARY, fontsize=9)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False,
               fontsize=8, labelcolor=TEXT_SECONDARY, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        f"Clade C subclades — distance from target centre to nearest {boundary} codon\n"
        "grey = matched random-placement null. Several subclades have n < 30: descriptive only.",
        color=TEXT_PRIMARY, fontsize=12,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.92))
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {out_path}")


# --------------------------------------------------------------------------

CAVEATS = """# Caveats for the per-clade IS1182 target analysis

These apply to `clade_summary.csv`, `ks_tests.csv` and the ECDF figures in this
directory. They are properties of the input data, not of the pipeline.

0. **`IS1182_all` is the union of `cladeA` + `cladeB` + `cladeC_all`.** It is
   plotted alongside them for reference, but it is a whole shown with its own
   parts — never treat it as a fourth independent group, and never test a clade
   against it. `IS4_all` is a different IS family and *is* disjoint from all of
   them; `ks_tests.csv` only ever compares the pairs listed in `DISJOINT_PAIRS`.

1. **Clade C subclades are nested inside `cladeC_all`.** 686 of the 875
   `cladeC_all` targets also appear in one of the seven subclade files, and 189
   `cladeC_all` targets have no subclade assignment at all. Never treat a
   subclade and `cladeC_all` as independent samples, and never sum them.
   Only `cladeA`, `cladeB` and `cladeC_all` are mutually disjoint — together
   they are exactly the 4121 targets in `IS1182_all`, which is therefore also
   not independent of them.

2. **One target is assigned to two subclades.** `CP181075_4693439-4693380`
   appears in both `cladeC_pcc4` and `cladeC_unknown5`.

3. **Several subclades are far too small for inference.** n = 10
   (`unknown5`), 15 (`unknown3`), 20 (`fba3`), 21 (`unknown4`), 30
   (`unknown1`). KS tests at these sizes are badly underpowered; read those
   panels descriptively and rely on the reported n, not the p-value.

4. **Targets are pseudoreplicated.** Multiple targets come from the same genome
   (up to 64 from one genome in cladeA), and targets in closely related genomes
   may be the same ancestral insertion counted many times. Two families in
   `ks_tests.csv` guard against this by using one value per accession (the
   per-genome median): `real_vs_null_genome_collapsed` (unpaired KS) and
   `real_vs_null_paired_by_genome` (Wilcoxon signed-rank). The paired test is
   the one that matches the design — the null is generated inside each genome,
   so real and null are paired by genome, and pairing also removes
   between-genome variation in gene density. **Prefer the paired result**; where
   it disagrees with the uncollapsed KS, the uncollapsed KS is inflated by
   pseudoreplication. A negative `median_diff_real_minus_null` means observed
   targets sit closer to the codon than random placements in the same genomes.

   `real_vs_null_sign_test_by_genome` reports the same pairs as a sign test —
   simply *how many* genomes show the effect. Read it alongside the Wilcoxon:
   these distances are bounded below by 0 and unbounded above, so a minority of
   genomes whose targets sit very far from any codon can dominate the
   magnitude-weighted signed ranks and cancel a consistent but modest shift in
   the majority. **`IS4_all` start/downstream is exactly this case**: 59% of its
   genomes sit closer than their own null (sign test, Holm p ~ 1e-15), while the
   Wilcoxon is not significant. Neither number is wrong. For stop/upstream the
   two tests agree, but IS4 is still far more heterogeneous between genomes
   (67% of genomes closer than their null) than IS1182, whose clades run at
   80-96%.

5. **The same locus is often present twice, once per orientation.** A forward
   entry (`X_100-159`) and a reverse entry (`X_159-100`) describe the same 60 bp
   window and yield distances 1 bp apart, so that insertion site is counted
   twice. `n_analysed_duplicate_windows` in `clade_summary.csv` reports how many
   analysed rows are the second member of such a pair (several hundred in the
   larger clades). They are counted, not merged, because whether the pair
   represents two independent insertions or one site recorded twice is a
   question about how the target set was built, not something this script can
   decide. If they are redundant, collapse on
   `accession` + `seq_start` + `seq_end` before drawing conclusions.

6. **Gene boundaries come from both `gene` and `CDS` features.** Duplicate
   gene/CDS pairs are harmless because only the minimum distance is used, but
   roughly 2% of `gene` features are tRNA/rRNA/tmRNA/ncRNA with no CDS and are
   credited with a "start/stop codon" they do not have. This was left unchanged
   so the clade numbers stay comparable with the earlier pooled runs.

7. **Distances are measured from the centre of the trimmed target**
   (position 31 of 60, in the target's own stranded orientation), mapped through
   the alignment. The random-placement null uses the same anchor offset, so the
   two are directly comparable. This differs from the earlier pooled runs in
   `results/IS1182_all_*`, which measured from whichever alignment endpoint sat
   closest to a gene — a min-of-two statistic that shifts distances downward.
   The two sets of numbers are not interchangeable.

8. **Coordinates in `coordinates_with_genes.csv` are 1-based inclusive.**
   Distances are computed 0-based internally against BioPython feature
   coordinates, so `up_dist`/`down_dist` = 0 means the anchor nucleotide is the
   boundary nucleotide. The earlier pooled runs carried a 1 bp offset in each
   direction; do not compare distances across the two versions at bp resolution.

9. **A handful of IS4 targets are not on bacterial genomes.** Four `IS4_all`
   headers carry RefSeq transcript/region accessions — `NG_032198` (human
   pseudogene), `NM_001126745` (*Xenopus* mRNA), `NM_105477` (*Arabidopsis*
   mRNA), `XM_073787841` (dolphin predicted mRNA). They are spurious for an IS
   target analysis. Earlier runs never showed them because the FASTA header
   parser silently dropped any accession containing an underscore; now that the
   parser is fixed they are processed and counted in
   `n_analysed_non_genomic`. They are 4 of 9956 records and move nothing, but
   they are left in rather than silently filtered — drop them upstream if the
   target set is regenerated.

10. **Genomes with no usable annotation are dropped, not counted as zero.**
   Rows with no flanking gene on either side are filtered out before any
   distance statistic; `n_unannotated` in `clade_summary.csv` reports how many.
   In practice these are GenBank records with no `gene`/`CDS` features at all,
   so nothing can be measured for them; re-annotating those genomes is the
   only way to bring them in.

11. **Targets inside a gene are counted.** Every annotated target contributes a
   distance, whether it sits inside a gene, overlaps a gene edge, or is
   intergenic (`n_inside`, `n_partial`, `n_intergenic`). Each distance runs to
   the nearest boundary on that side — the gene the target sits in, or a
   neighbouring gene if its boundary is nearer (`anchor_in_gene` vs
   `up_gene`/`down_gene` in `coordinates_with_genes.csv`). The random null keeps
   inside-gene placements too, so real and null still match. The previous
   analysis in `results/IS1182_clades` dropped every target fully inside a gene
   (and every null placement inside a gene), which removed the sites furthest
   from a codon and pulled both ECDFs towards zero; its numbers are not
   comparable with these. `comparison_vs_v1.csv` tabulates the change, and
   `ecdf_by_location_*.png` shows each location class separately — the
   all-target curve is a mixture of them, so a shift in it can come from the mix
   as well as from the distances.

12. **Origin-spanning genes and circular molecules.** A gene crossing the origin
   of a circular molecule is handled as its real segments, with codons taken
   from its first and last parts. Previously such a gene spanned the entire
   chromosome, so every target in that genome was classified inside a gene and
   dropped — true intergenic targets included — and its codons sat at the
   genome ends. On circular records, distances now also wrap the origin.
   `comparison_location_v1_v2.csv` shows how the old `inside` rows reclassify.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results_dir", required=True,
                    help="Directory holding the <clade>_<start|stop> run directories.")
    ap.add_argument("--fasta_dir", required=True,
                    help="Directory holding the per-clade target FASTA files.")
    ap.add_argument("--xmax", type=float, nargs="+", default=[400.0, 2000.0],
                    help="x-axis limit(s) for the ECDF figures (bp); one set of figures "
                         "per value, suffixed _x<N>. Default 400 2000.")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    fasta_dir = Path(args.fasta_dir)

    # ---- discover runs -------------------------------------------------
    data: dict[tuple[str, str], dict] = {}
    for d in sorted(results_dir.iterdir()):
        m = re.fullmatch(r"(.+)_(start|stop)", d.name) if d.is_dir() else None
        if not m:
            continue
        clade, boundary = m.group(1), m.group(2)
        real, null = load_run(d)
        if real is None:
            print(f"WARNING: {d.name} has no coordinates_with_genes.csv — skipped")
            continue
        data[(clade, boundary)] = {"dir": d, "raw": real,
                                   "real": analysed_subset(real), "null": null}
    if not data:
        raise SystemExit(f"No <clade>_<boundary> result directories found in {results_dir}")
    print(f"Loaded {len(data)} runs: "
          f"{sorted({c for c, _ in data})}")

    # ---- summary table -------------------------------------------------
    rows = []
    for (clade, boundary), e in sorted(data.items()):
        raw, real, null = e["raw"], e["real"], e["null"]
        n_input = count_fasta_records(fasta_dir / FASTA_FOR.get(clade, "___missing___"))
        ok = raw[raw["error"].isna()]
        n_unannotated = int((ok["up_dist"].isna() & ok["down_dist"].isna()).sum())
        # Location of analysed rows only, so inside + partial + intergenic == n_analysed.
        loc = real["location"].value_counts()
        for direction in DIRECTIONS:
            r = describe(real[f"{direction}_dist"])
            base = dict(
                clade=clade, boundary_type=boundary, direction=direction,
                n_fasta_records=n_input,
                n_rows=len(raw),
                n_parse_dropped=(n_input - len(raw)) if n_input >= 0 else np.nan,
                n_error=int(raw["error"].notna().sum()),
                n_inside=int(loc.get("inside", 0)),
                n_partial=int(loc.get("partial", 0)),
                n_intergenic=int(loc.get("intergenic", 0)),
                n_unannotated=n_unannotated,
                n_analysed=len(real),
                n_anchor_in_gene=int(real["anchor_in_gene"].notna().sum())
                if "anchor_in_gene" in real.columns else np.nan,
                n_genomes=int(real["accession"].nunique()),
                # Forward and reverse entries at the same locus normalise to the
                # same 60 bp window, so they contribute two near-identical
                # distances for one insertion site. Counted, not silently merged.
                n_analysed_duplicate_windows=int(
                    len(real) - real[["accession", "seq_start", "seq_end"]]
                    .drop_duplicates().shape[0]
                ),
                # RefSeq transcript/region accessions are not bacterial genomes;
                # any hit on one is spurious for an IS target analysis.
                n_analysed_non_genomic=int(
                    real["accession"].astype(str)
                    .str.match(r"^(NM_|XM_|XR_|NR_|NG_)").sum()
                ),
                n_midpoint_fallback=int((raw["anchor_source"] == "midpoint_fallback").sum())
                if "anchor_source" in raw.columns else np.nan,
            )
            base.update({f"real_{k}": v for k, v in r.items()})
            if null is not None:
                rn = describe(null[f"{direction}_dist"])
                base.update({f"null_{k}": v for k, v in rn.items()})
                for cut in (50, 100, 200):
                    denom = rn[f"frac_le{cut}"]
                    base[f"enrichment_le{cut}"] = (
                        r[f"frac_le{cut}"] / denom if denom and denom > 0 else np.nan
                    )
            rows.append(base)

    summary = pd.DataFrame(rows)
    summary_path = results_dir / "clade_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"wrote {summary_path}  ({len(summary)} rows)")

    # ---- KS tests ------------------------------------------------------
    tests = []

    # (a) real vs its own null, per clade / boundary / direction
    for (clade, boundary), e in sorted(data.items()):
        if e["null"] is None:
            continue
        for direction in DIRECTIONS:
            col = f"{direction}_dist"
            D, p, n1, n2 = ks(e["real"][col], e["null"][col])
            tests.append(dict(family="real_vs_null", clade=clade, other="null",
                              boundary_type=boundary, direction=direction,
                              n1=n1, n2=n2, statistic=D, p=p))
            # genome-collapsed: one value per accession on both sides, so that a
            # genome contributing 64 targets counts once rather than 64 times.
            rc = e["real"].groupby("accession")[col].median()
            nc = e["null"].groupby("accession")[col].median()
            D, p, n1, n2 = ks(rc, nc)
            tests.append(dict(family="real_vs_null_genome_collapsed", clade=clade,
                              other="null", boundary_type=boundary, direction=direction,
                              n1=n1, n2=n2, statistic=D, p=p))
            # ...and the same genomes paired, which is what the design actually is
            W, p, npairs, mdiff = paired_wilcoxon(rc, nc)
            tests.append(dict(family="real_vs_null_paired_by_genome", clade=clade,
                              other="null", boundary_type=boundary, direction=direction,
                              n1=npairs, n2=npairs, statistic=W, p=p,
                              median_diff_real_minus_null=mdiff))
            # Sign test on the same pairs: robust to the heavy right tail that
            # can cancel the Wilcoxon signed ranks (see paired_sign_test).
            frac, p, nsign = paired_sign_test(rc, nc)
            tests.append(dict(family="real_vs_null_sign_test_by_genome", clade=clade,
                              other="null", boundary_type=boundary, direction=direction,
                              n1=nsign, n2=nsign, statistic=frac, p=p,
                              median_diff_real_minus_null=mdiff,
                              frac_genomes_real_below_null=frac))

    # (b) pairwise between sets that are genuinely disjoint (see DISJOINT_PAIRS)
    for boundary in BOUNDARIES:
        for direction in DIRECTIONS:
            col = f"{direction}_dist"
            for a, b in DISJOINT_PAIRS:
                ea, eb = data.get((a, boundary)), data.get((b, boundary))
                if ea is None or eb is None:
                    continue
                D, p, n1, n2 = ks(ea["real"][col], eb["real"][col])
                tests.append(dict(family="set_vs_set", clade=a, other=b,
                                  boundary_type=boundary, direction=direction,
                                  n1=n1, n2=n2, statistic=D, p=p))
                D, p, n1, n2 = ks(ea["real"].groupby("accession")[col].median(),
                                  eb["real"].groupby("accession")[col].median())
                tests.append(dict(family="set_vs_set_genome_collapsed",
                                  clade=a, other=b, boundary_type=boundary,
                                  direction=direction, n1=n1, n2=n2, statistic=D, p=p))

    ks_df = pd.DataFrame(tests)
    if not ks_df.empty:
        # Holm correction applied within each test family, not across all of them.
        ks_df["p_holm"] = np.nan
        for fam, grp in ks_df.groupby("family"):
            ks_df.loc[grp.index, "p_holm"] = holm(list(grp["p"]))
        ks_df = ks_df.sort_values(["family", "boundary_type", "direction", "clade"])
    ks_path = results_dir / "ks_tests.csv"
    ks_df.to_csv(ks_path, index=False)
    print(f"wrote {ks_path}  ({len(ks_df)} tests)")

    write_prism_cdf(data, results_dir / "prism_cdf_clades.csv")

    # ---- figures -------------------------------------------------------
    for xmax in args.xmax:
        sfx = f"_x{int(xmax)}"
        plot_series_grid(
            data, FAMILIES, xmax, results_dir / f"ecdf_families{sfx}.png",
            "IS4 vs IS1182 target sites: distance from target centre to nearest flanking codon",
            ncol_legend=2,
        )
        plot_series_grid(
            data, MAIN_SERIES, xmax, results_dir / f"ecdf_main_clades{sfx}.png",
            "IS4 and IS1182 target sites, whole families and IS1182 clades",
            subtitle="IS1182_all is the union of cladeA + cladeB + cladeC_all — a whole "
                     "and its parts, not independent series",
        )
        for fam in FAMILIES:
            if any((fam, b) in data for b in BOUNDARIES):
                plot_by_location(data, fam, xmax,
                                 results_dir / f"ecdf_by_location_{fam}{sfx}.png")
        for boundary in BOUNDARIES:
            if any((c, boundary) in data for c in SUBCLADES):
                plot_subclades(data, boundary, xmax,
                               results_dir / f"ecdf_cladeC_subclades_{boundary}{sfx}.png")

    (results_dir / "CAVEATS.md").write_text(CAVEATS)
    print(f"wrote {results_dir / 'CAVEATS.md'}")

    # ---- console digest ------------------------------------------------
    print("\n=== headline: median distance from target centre (bp), real vs null ===")
    cols = ["clade", "boundary_type", "direction", "n_analysed", "n_inside",
            "real_n", "real_median", "null_median", "enrichment_le100"]
    have = [c for c in cols if c in summary.columns]
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summary[have].to_string(index=False))

    if not ks_df.empty:
        print("\n=== real vs matched null, paired by genome (Holm-adjusted) ===")
        print("negative median_diff = observed targets sit CLOSER to the codon "
              "than random placements in the same genomes;")
        print("frac_below = fraction of genomes where that holds (sign test).")
        w = ks_df[ks_df["family"] == "real_vs_null_paired_by_genome"].set_index(
            ["clade", "boundary_type", "direction"])
        s = ks_df[ks_df["family"] == "real_vs_null_sign_test_by_genome"].set_index(
            ["clade", "boundary_type", "direction"])
        merged = w[["n1", "median_diff_real_minus_null", "p_holm"]].join(
            s[["frac_genomes_real_below_null", "p_holm"]],
            lsuffix="_wilcoxon", rsuffix="_sign").reset_index()
        with pd.option_context("display.width", 220, "display.max_rows", 400):
            print(merged.round(4).to_string(index=False))

        print("\n=== real vs matched null (KS, Holm-adjusted within family) ===")
        sel = ks_df[ks_df["family"].str.startswith("real_vs_null")]
        with pd.option_context("display.width", 200, "display.max_rows", 400):
            print(sel[["family", "clade", "boundary_type", "direction",
                       "n1", "n2", "statistic", "p_holm"]].to_string(index=False))

    print("\nSee CAVEATS.md before interpreting any of the above.")


if __name__ == "__main__":
    main()
