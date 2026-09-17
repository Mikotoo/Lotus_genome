#!/usr/bin/env python3
"""
Gifu telomere analysis (motif calling, coverage, spanning reads, plots).

Use subcommands: find-tel, coverage, alignments, plot, all

Steps:
  1. Scan 7-bp telomere units (TTTAGGG / CCCTAAA, <=1 mismatch) on the plus strand.
  2. Compute 100-bp window density as sequence coverage (merged unit spans / bin width).
  3. Call telomere bins from each chromosome end inward; bridge up to 5 consecutive
     low-density bins when high-density bins continue on both sides.
  4. Refine coordinates from the outermost and innermost telomere bins: walk
     consecutive units from each bin edge toward the chromosome end; final span
     is the union of both walks (gaps in the middle are ignored).
  5. Write coordinates TSV; optional PDF figures if matplotlib is available.

Units are deduplicated by (chrom, start) before density and coordinate calling.
"""

from __future__ import annotations

import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
import argparse
import bisect
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

MOTIFS = ("CCCTAAA", "TTTAGGG")
BASES = "ATGC"
UNIT_LEN = 7
DEFAULT_BIN = 100
DEFAULT_DENSITY = 0.5
DEFAULT_MAX_UNIT_GAP = 14
DEFAULT_MAX_DIP_BINS = 5
DEFAULT_PLOT_FLANK = 10_000
MAX_DENSITY = 1.0
MOTIF_COLORS = {"CCCTAAA": "#2E86AB", "TTTAGGG": "#C73E1D"}


@dataclass
class TelUnit:
    chrom: str
    start: int
    end: int
    motif: str


@dataclass
class BinRow:
    chrom: str
    bin_start: int
    bin_end: int
    count: int
    density: float


@dataclass
class TelomereCall:
    chrom: str
    end: str
    bin_start: int
    bin_end: int
    tel_start: int
    tel_end: int
    length_bp: int
    dominant_motif: str
    n_units: int
    mean_bin_density: float
    max_bin_density: float


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    chrom: Optional[str] = None
    chunks: List[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    yield chrom, "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if chrom is not None:
            yield chrom, "".join(chunks).upper()


def scan_telomere_units(fasta_path: Path) -> List[TelUnit]:
    units: List[TelUnit] = []
    for chrom, seq in parse_fasta(fasta_path):
        n = len(seq)
        if n < UNIT_LEN:
            continue
        for motif in MOTIFS:
            for i in range(n - UNIT_LEN + 1):
                for j in range(UNIT_LEN):
                    for base in BASES:
                        variant = motif[:j] + base + motif[j + 1 :]
                        if seq[i : i + UNIT_LEN] == variant:
                            units.append(TelUnit(chrom, i, i + UNIT_LEN, variant))
    return dedupe_units(units)


def dedupe_units(units: Sequence[TelUnit]) -> List[TelUnit]:
    """One record per (chrom, start); prefer exact TTTAGGG / CCCTAAA motifs."""
    best: Dict[Tuple[str, int], TelUnit] = {}
    for u in units:
        key = (u.chrom, u.start)
        cur = best.get(key)
        if cur is None:
            best[key] = u
            continue
        if cur.motif in MOTIFS and u.motif not in MOTIFS:
            continue
        if u.motif in MOTIFS:
            best[key] = u
    return sorted(best.values(), key=lambda u: (u.chrom, u.start))


def write_unit_bed(units: Sequence[TelUnit], bed_path: Path) -> None:
    with bed_path.open("w") as fh:
        for u in units:
            fh.write(f"{u.chrom}\t{u.start}\t{u.end}\t{u.motif}\n")


def load_units_bed(path: Path) -> List[TelUnit]:
    units: List[TelUnit] = []
    with path.open() as fh:
        for line in fh:
            c, s, e, m = line.rstrip("\n").split("\t")[:4]
            start, end = int(s), int(e)
            if end - start == 6:
                end = start + UNIT_LEN
            units.append(TelUnit(c, start, end, m))
    return dedupe_units(units)


def load_chrom_lengths(fai_path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with fai_path.open() as fh:
        for line in fh:
            parts = line.split()
            lengths[parts[0]] = int(parts[1])
    return lengths


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for lo, hi in intervals[1:]:
        prev_lo, prev_hi = merged[-1]
        if lo <= prev_hi:
            merged[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged.append((lo, hi))
    return merged


def _units_overlapping_bin(
    chr_units: Sequence[TelUnit], starts: Sequence[int], bin_lo: int, bin_hi: int
) -> List[TelUnit]:
    """Units whose 7-bp span overlaps [bin_lo, bin_hi)."""
    lo = bisect.bisect_left(starts, bin_lo - (UNIT_LEN - 1))
    hi = bisect.bisect_left(starts, bin_hi)
    return [u for u in chr_units[lo:hi] if u.start < bin_hi and u.end > bin_lo]


def _coverage_from_units(units: Sequence[TelUnit], bin_lo: int, bin_hi: int) -> float:
    """Fraction of bin width covered by merged telomere unit intervals (0–1)."""
    width = bin_hi - bin_lo
    if width <= 0:
        return 0.0
    intervals = [
        (max(bin_lo, u.start), min(bin_hi, u.end))
        for u in units
        if u.start < bin_hi and u.end > bin_lo
    ]
    if not intervals:
        return 0.0
    covered_bp = sum(hi - lo for lo, hi in _merge_intervals(intervals))
    return min(MAX_DENSITY, covered_bp / width)


def _motif_coverage_in_bin(
    units: Sequence[TelUnit], bin_lo: int, bin_hi: int
) -> Tuple[float, float]:
    """Merged-sequence coverage for exact CCCTAAA and TTTAGGG units separately."""
    width = bin_hi - bin_lo
    if width <= 0:
        return 0.0, 0.0
    ccct_iv: List[Tuple[int, int]] = []
    ttta_iv: List[Tuple[int, int]] = []
    for u in units:
        if u.start < bin_hi and u.end > bin_lo:
            seg = (max(bin_lo, u.start), min(bin_hi, u.end))
            if u.motif == "CCCTAAA":
                ccct_iv.append(seg)
            elif u.motif == "TTTAGGG":
                ttta_iv.append(seg)
    ccct_cov = (
        sum(hi - lo for lo, hi in _merge_intervals(ccct_iv)) / width if ccct_iv else 0.0
    )
    ttta_cov = (
        sum(hi - lo for lo, hi in _merge_intervals(ttta_iv)) / width if ttta_iv else 0.0
    )
    return min(MAX_DENSITY, ccct_cov), min(MAX_DENSITY, ttta_cov)


def compute_bin_density(
    units: Sequence[TelUnit],
    chrom_lengths: Dict[str, int],
    bin_size: int,
) -> List[BinRow]:
    by_chr: Dict[str, List[TelUnit]] = {}
    for u in units:
        by_chr.setdefault(u.chrom, []).append(u)
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda u: u.start)

    rows: List[BinRow] = []
    for chrom, length in chrom_lengths.items():
        chr_units = by_chr.get(chrom, [])
        starts = [u.start for u in chr_units]
        n_bins = (length + bin_size - 1) // bin_size
        for b in range(n_bins):
            bs = b * bin_size
            be = min(bs + bin_size, length)
            subset = _units_overlapping_bin(chr_units, starts, bs, be)
            count = len(subset)
            density = _coverage_from_units(subset, bs, be)
            rows.append(BinRow(chrom, bs, be, count, density))
    return rows


def _motif_density_in_bin(
    units: Sequence[TelUnit], bin_lo: int, bin_hi: int, bin_size: int
) -> Tuple[float, float]:
    del bin_size  # width uses actual bin span
    subset = [u for u in units if u.start < bin_hi and u.end > bin_lo]
    return _motif_coverage_in_bin(subset, bin_lo, bin_hi)


def _plot_panel_bins(
    units: Sequence[TelUnit],
    chrom: str,
    chrom_len: int,
    end: str,
    flank_bp: int,
    bin_size: int,
) -> List[Tuple[float, float, float]]:
    """
    Bins within flank_bp of one chromosome end.

    Returns (x_from_end, ccct_density, tttaggg_density) per bin.
    x_from_end is 0 at the telomere end and increases inward.
    """
    chr_units = [u for u in units if u.chrom == chrom]
    if end == "5p":
        lo, hi = 0, min(flank_bp, chrom_len)
    else:
        lo, hi = max(0, chrom_len - flank_bp), chrom_len

    out: List[Tuple[float, float, float]] = []
    b = (lo // bin_size) * bin_size
    while b < hi:
        be = min(b + bin_size, chrom_len)
        if be <= b:
            break
        center = (b + be) / 2
        if end == "5p":
            x_from_end = center
        else:
            x_from_end = chrom_len - center
        ccct, ttta = _motif_density_in_bin(chr_units, b, be, bin_size)
        out.append((x_from_end, ccct, ttta))
        b += bin_size
    return out


def _bins_for_chrom(rows: Sequence[BinRow], chrom: str) -> List[BinRow]:
    return sorted(
        [r for r in rows if r.chrom == chrom],
        key=lambda r: r.bin_start,
    )


def _contiguous_tel_bins(
    bins: Sequence[BinRow],
    from_5p: bool,
    density_threshold: float,
    max_dip_bins: int = DEFAULT_MAX_DIP_BINS,
) -> List[BinRow]:
    """
    Extend telomere bins from the chromosome end inward.

    Leading bins at the physical chromosome end with density <= threshold are
    skipped (common for partial terminal bins). Extension starts at the first
    high-density bin encountered when moving inward.

    A run of bins with density <= threshold may be bridged (included) when it is
    at most max_dip_bins long and is flanked by high-density bins on both sides
    (high-density from the chromosome end before the dip, and recovery after).
    """
    if not bins:
        return []
    ordered = list(bins) if from_5p else list(reversed(bins))
    selected: List[BinRow] = []
    i = 0
    n = len(ordered)
    while i < n:
        row = ordered[i]
        if row.density > density_threshold:
            selected.append(row)
            i += 1
            continue
        if not selected:
            i += 1
            continue
        dip_start = i
        while i < n and ordered[i].density <= density_threshold:
            i += 1
        dip_len = i - dip_start
        if dip_len > max_dip_bins:
            break
        if i >= n or ordered[i].density <= density_threshold:
            break
        selected.extend(ordered[dip_start:i])
    return sorted(selected, key=lambda r: r.bin_start)


def _units_in_bin(
    units: Sequence[TelUnit], bin_lo: int, bin_hi: int
) -> List[TelUnit]:
    return [u for u in units if bin_lo <= u.start < bin_hi]


def _walk_from_seed(
    region_units: Sequence[TelUnit],
    seed: TelUnit,
    max_unit_gap: int,
    toward_increasing: bool,
) -> List[TelUnit]:
    """Extend from seed along consecutive units in one direction."""
    idx = next(i for i, u in enumerate(region_units) if u.start == seed.start)
    cluster = [seed]
    if toward_increasing:
        for i in range(idx + 1, len(region_units)):
            if region_units[i].start - cluster[-1].start <= max_unit_gap:
                cluster.append(region_units[i])
            else:
                break
    else:
        for i in range(idx - 1, -1, -1):
            if cluster[0].start - region_units[i].start <= max_unit_gap:
                cluster.insert(0, region_units[i])
            else:
                break
    return cluster


def _refine_from_outer_inner_bins(
    units: Sequence[TelUnit],
    from_5p: bool,
    region_lo: int,
    region_hi: int,
    tel_bins: Sequence[BinRow],
    max_unit_gap: int,
) -> List[TelUnit]:
    """
    From the outermost and innermost telomere bins, walk consecutive units at
    each bin edge toward the chromosome end. Return the union of both walks;
    interior gaps between the two walks are ignored for the final span.
    """
    region_units = sorted(
        [u for u in units if region_lo <= u.start < region_hi],
        key=lambda u: u.start,
    )
    if not region_units or not tel_bins:
        return []

    if from_5p:
        outer_bin, inner_bin = tel_bins[0], tel_bins[-1]
        outer_units = _units_in_bin(
            region_units, outer_bin.bin_start, outer_bin.bin_end
        )
        inner_units = _units_in_bin(
            region_units, inner_bin.bin_start, inner_bin.bin_end
        )
        if not outer_units or not inner_units:
            return []
        outer_seed = min(outer_units, key=lambda u: u.start)
        inner_seed = max(inner_units, key=lambda u: u.start)
        outer_walk = _walk_from_seed(
            region_units, outer_seed, max_unit_gap, toward_increasing=True
        )
        inner_walk = _walk_from_seed(
            region_units, inner_seed, max_unit_gap, toward_increasing=False
        )
    else:
        outer_bin, inner_bin = tel_bins[-1], tel_bins[0]
        outer_units = _units_in_bin(
            region_units, outer_bin.bin_start, outer_bin.bin_end
        )
        inner_units = _units_in_bin(
            region_units, inner_bin.bin_start, inner_bin.bin_end
        )
        if not outer_units or not inner_units:
            return []
        outer_seed = max(outer_units, key=lambda u: u.start)
        inner_seed = min(inner_units, key=lambda u: u.start)
        outer_walk = _walk_from_seed(
            region_units, outer_seed, max_unit_gap, toward_increasing=True
        )
        inner_walk = _walk_from_seed(
            region_units, inner_seed, max_unit_gap, toward_increasing=True
        )

    seen = {u.start: u for u in outer_walk + inner_walk}
    return sorted(seen.values(), key=lambda u: u.start)


def call_telomere_at_end(
    chrom: str,
    chrom_len: int,
    bin_rows: Sequence[BinRow],
    units: Sequence[TelUnit],
    end: str,
    density_threshold: float,
    max_unit_gap: int,
    max_dip_bins: int,
) -> Optional[TelomereCall]:
    from_5p = end == "5p"
    tel_bins = _contiguous_tel_bins(
        bin_rows, from_5p, density_threshold, max_dip_bins
    )
    if not tel_bins:
        return None

    bin_start = min(b.bin_start for b in tel_bins)
    bin_end = max(b.bin_end for b in tel_bins)

    chr_units = [u for u in units if u.chrom == chrom]
    if from_5p:
        region_lo, region_hi = 0, bin_end
    else:
        region_lo, region_hi = bin_start, chrom_len

    walk_units = _refine_from_outer_inner_bins(
        chr_units,
        from_5p,
        region_lo,
        region_hi,
        tel_bins,
        max_unit_gap,
    )
    if not walk_units:
        return None

    tel_start = min(u.start for u in walk_units)
    tel_end = max(u.end for u in walk_units)
    span_units = [
        u
        for u in chr_units
        if region_lo <= u.start < region_hi and tel_start <= u.start < tel_end
    ]
    dominant = Counter(u.motif for u in span_units).most_common(1)[0][0]
    densities = [min(MAX_DENSITY, b.density) for b in tel_bins]

    return TelomereCall(
        chrom=chrom,
        end=end,
        bin_start=bin_start,
        bin_end=bin_end,
        tel_start=tel_start,
        tel_end=tel_end,
        length_bp=tel_end - tel_start,
        dominant_motif=dominant,
        n_units=len(span_units),
        mean_bin_density=sum(densities) / len(densities),
        max_bin_density=max(densities),
    )


def call_all_telomeres(
    bin_rows: Sequence[BinRow],
    units: Sequence[TelUnit],
    chrom_lengths: Dict[str, int],
    density_threshold: float,
    max_unit_gap: int,
    max_dip_bins: int,
) -> List[TelomereCall]:
    calls: List[TelomereCall] = []
    for chrom, chrom_len in chrom_lengths.items():
        sub = _bins_for_chrom(bin_rows, chrom)
        for end in ("5p", "3p"):
            call = call_telomere_at_end(
                chrom,
                chrom_len,
                sub,
                units,
                end,
                density_threshold,
                max_unit_gap,
                max_dip_bins,
            )
            if call:
                calls.append(call)
    return calls


def write_density_table(bin_rows: Sequence[BinRow], path: Path) -> None:
    with path.open("w") as fh:
        for r in bin_rows:
            wid = f"{r.chrom}_{r.bin_start}_{r.bin_end}"
            fh.write(f"{wid}\t{r.density:.4g}\n")


def write_density_detail(bin_rows: Sequence[BinRow], path: Path) -> None:
    fields = ["chrom", "bin_start", "bin_end", "count", "density"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in bin_rows:
            w.writerow(
                {
                    "chrom": r.chrom,
                    "bin_start": r.bin_start,
                    "bin_end": r.bin_end,
                    "count": r.count,
                    "density": f"{r.density:.4g}",
                }
            )


def write_coords(calls: Sequence[TelomereCall], path: Path) -> None:
    fields = [
        "chromosome",
        "end",
        "bin_region_start",
        "bin_region_end",
        "telomere_start",
        "telomere_end",
        "length_bp",
        "dominant_motif",
        "n_units",
        "mean_bin_density",
        "max_bin_density",
    ]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for c in calls:
            w.writerow(
                {
                    "chromosome": c.chrom,
                    "end": c.end,
                    "bin_region_start": c.bin_start,
                    "bin_region_end": c.bin_end,
                    "telomere_start": c.tel_start,
                    "telomere_end": c.tel_end,
                    "length_bp": c.length_bp,
                    "dominant_motif": c.dominant_motif,
                    "n_units": c.n_units,
                    "mean_bin_density": f"{c.mean_bin_density:.4f}",
                    "max_bin_density": f"{c.max_bin_density:.4f}",
                }
            )


def _natural_key(chrom: str):
    m = re.search(r"(\d+)$", chrom)
    return int(m.group(1)) if m else chrom


def _configure_plot_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _format_bp(pos: int) -> str:
    return f"{pos:,}"


def _tel_span_plot_coords(
    call: Optional[TelomereCall],
    chrom_len: int,
    end: str,
    flank_bp: int,
) -> Optional[Tuple[float, float]]:
    """Telomere span as distance (bp) from the chromosome end, within the plot window."""
    if call is None:
        return None
    if end == "5p":
        x0, x1 = float(call.tel_start), float(call.tel_end)
    else:
        x0 = chrom_len - call.tel_end
        x1 = chrom_len - call.tel_start
    x0 = max(0.0, x0)
    x1 = min(float(flank_bp), x1)
    if x1 <= x0:
        return None
    return x0, x1


def _draw_end_column(
    fig,
    col_gs,
    chrom: str,
    chrom_len: int,
    end: str,
    call: Optional[TelomereCall],
    panel_bins: Sequence[Tuple[float, float, float]],
    flank_bp: int,
    bin_size: int,
    density_threshold: float,
) -> None:
    """
    One chromosome end: shared x-axis; top = telomere marker + x ticks;
    bottom = stacked motif density (no bottom x ticks).
    """
    inner = col_gs.subgridspec(2, 1, height_ratios=[0.70, 3.0], hspace=0.06)
    ax_top = fig.add_subplot(inner[0])
    ax_bot = fig.add_subplot(inner[1], sharex=ax_top)

    ax_top.set_xlim(0, flank_bp)
    if end == "3p":
        ax_top.invert_xaxis()

    ax_top.set_ylim(0, 1)
    ax_top.set_yticks([])
    span = _tel_span_plot_coords(call, chrom_len, end, flank_bp)
    if span is not None and call is not None:
        x0, x1 = span
        ax_top.plot(
            [x0, x1],
            [0.38, 0.38],
            color="#C73E1D",
            solid_capstyle="butt",
            linewidth=3.5,
            clip_on=False,
            zorder=3,
        )
        ax_top.text(
            (x0 + x1) / 2,
            0.58,
            f"{_format_bp(call.tel_start)} – {_format_bp(call.tel_end)} bp",
            ha="center",
            va="bottom",
            fontsize=7,
            fontfamily="Liberation Sans",
            zorder=4,
        )
    elif call is None:
        ax_top.text(
            flank_bp / 2,
            0.45,
            "No telomere call",
            ha="center",
            va="center",
            fontsize=7,
            fontfamily="Liberation Sans",
            color="#888888",
        )

    ax_top.tick_params(
        axis="x",
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        length=4,
        width=0.8,
        labelsize=7,
        direction="in",
    )
    for lbl in ax_top.get_xticklabels():
        lbl.set_fontfamily("Liberation Sans")
    ax_top.spines["top"].set_visible(True)
    ax_top.spines["top"].set_linewidth(0.6)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.set_title(f"{chrom} {end}", fontsize=8, fontfamily="Liberation Sans", pad=3)

    if not panel_bins:
        ax_bot.text(
            0.5,
            0.5,
            "No bins in window",
            transform=ax_bot.transAxes,
            ha="center",
            fontfamily="Liberation Sans",
        )
    else:
        xs = [b[0] for b in panel_bins]
        ccct = [b[1] for b in panel_bins]
        ttta = [b[2] for b in panel_bins]
        width = bin_size * 0.92

        ax_bot.bar(
            xs,
            ccct,
            width=width,
            label="CCCTAAA",
            color=MOTIF_COLORS["CCCTAAA"],
            edgecolor="none",
            align="center",
        )
        ax_bot.bar(
            xs,
            ttta,
            width=width,
            bottom=ccct,
            label="TTTAGGG",
            color=MOTIF_COLORS["TTTAGGG"],
            edgecolor="none",
            align="center",
        )
        ax_bot.axhline(
            density_threshold, color="#888888", ls="--", lw=0.9, dashes=(4, 3)
        )

    ax_bot.set_ylim(0, 1.05)
    ax_bot.tick_params(
        axis="x",
        labelbottom=False,
        bottom=False,
        top=False,
    )
    ax_bot.set_xlabel(
        "Distance from telomere end (bp)", fontsize=7, fontfamily="Liberation Sans", labelpad=3
    )
    ax_bot.set_ylabel("Unit density", fontsize=7, fontfamily="Liberation Sans")
    ax_bot.tick_params(axis="y", labelsize=7)
    for lbl in ax_bot.get_yticklabels():
        lbl.set_fontfamily("Liberation Sans")
    legend_loc = "upper left" if end == "3p" else "upper right"
    ax_bot.legend(
        loc=legend_loc,
        frameon=False,
        fontsize=7,
        prop={"family": "Liberation Sans"},
    )
    ax_bot.grid(False)
    for spine in ("top", "right"):
        ax_bot.spines[spine].set_visible(False)


def try_plot(
    units: Sequence[TelUnit],
    calls: Sequence[TelomereCall],
    chrom_lengths: Dict[str, int],
    outdir: Path,
    prefix: str,
    flank_bp: int,
    density_threshold: float,
    bin_size: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print(
            "matplotlib not found — skipping plots (coords TSV still written).",
            file=sys.stderr,
        )
        return

    _configure_plot_style()

    plot_dir = outdir / "plots"
    plot_dir.mkdir(exist_ok=True)
    calls_by = {(c.chrom, c.end): c for c in calls}

    chroms = sorted(chrom_lengths, key=_natural_key)
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    outer_gs = GridSpec(3, 2, figure=fig, hspace=0.25, wspace=0.35,
                        left=0.08, right=0.96, top=0.97, bottom=0.04)

    for idx, chrom in enumerate(chroms):
        chrom_len = chrom_lengths[chrom]
        row, col = idx // 2, idx % 2
        inner_gs = outer_gs[row, col].subgridspec(1, 2, wspace=0.18)

        for end_col, end in enumerate(("5p", "3p")):
            call = calls_by.get((chrom, end))
            panel_bins = _plot_panel_bins(
                units, chrom, chrom_len, end, flank_bp, bin_size
            )
            _draw_end_column(
                fig,
                inner_gs[end_col],
                chrom,
                chrom_len,
                end,
                call,
                panel_bins,
                flank_bp,
                bin_size,
                density_threshold,
            )

    fig.suptitle(
        f"Telomere ends (0–{_format_bp(flank_bp)} bp from each end) — {prefix}",
        fontsize=10,
        fontfamily="Liberation Sans",
        y=0.995,
    )
    fig.savefig(plot_dir / f"{prefix}_all_ends.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

    plot_telomere_lengths(calls, plot_dir, prefix)

    print(f"Figures saved under {plot_dir}/")


def plot_telomere_lengths(
    calls: Sequence[TelomereCall],
    plot_dir: Path,
    prefix: str,
) -> None:
    """Summary bar chart of telomere tract length per chromosome end."""
    import matplotlib.pyplot as plt

    if not calls:
        return

    plot_dir.mkdir(parents=True, exist_ok=True)
    _configure_plot_style()

    fig, ax = plt.subplots(figsize=(10, 5))
    labels = [f"{c.chrom}_{c.end}" for c in calls]
    lengths = [c.length_bp for c in calls]
    motifs = [c.dominant_motif for c in calls]
    x = list(range(len(calls)))
    bar_colors = [MOTIF_COLORS.get(m, "#888888") for m in motifs]
    ax.bar(x, lengths, color=bar_colors, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontfamily="Liberation Sans")
    ax.set_ylabel("Telomere tract length (bp)", fontfamily="Liberation Sans")
    ax.set_title("Telomere tract length per chromosome end", fontfamily="Liberation Sans")
    ax.grid(False)
    for i, (length, motif) in enumerate(zip(lengths, motifs)):
        ax.text(
            i,
            length + max(lengths) * 0.02,
            motif,
            ha="center",
            va="bottom",
            fontsize=8,
            fontfamily="Liberation Sans",
        )
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily("Liberation Sans")
    fig.tight_layout()
    fig.savefig(plot_dir / f"{prefix}_telomere_lengths.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {plot_dir}/{prefix}_telomere_lengths.pdf")


def resolve_fai(fasta: Path, fai_arg: Optional[str]) -> Path:
    """Resolve FASTA index path (follow symlinks)."""
    fasta = fasta.resolve()
    if fai_arg:
        fai = Path(fai_arg).resolve()
    else:
        fai = Path(str(fasta) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(
            f"FAI not found: {fai}\n"
            f"  Provide --fai or run: samtools faidx {fasta}"
        )
    return fai


def run_pipeline(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = Path(args.fasta).resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"FASTA not found: {fasta}")
    fai = resolve_fai(fasta, args.fai)
    prefix = args.prefix

    print(f"FASTA : {fasta}")
    print(f"FAI   : {fai}")
    chrom_lengths = load_chrom_lengths(fai)

    unit_path = None
    if args.reuse_units:
        unit_path = Path(args.reuse_units)
        if not unit_path.is_file():
            print(
                f"NOTE: unit bed not found ({unit_path}); "
                "scanning telomere units ..."
            )
            unit_path = None

    if unit_path is not None:
        print(f"Loading units from {unit_path} ...")
        with unit_path.open() as fh:
            raw_rows = sum(1 for line in fh if line.strip() and not line.startswith("#"))
        units = load_units_bed(unit_path)
        if len(units) < raw_rows:
            print(f"Deduplicating bed: {raw_rows} -> {len(units)} rows")
            write_unit_bed(units, unit_path)
    else:
        print(f"Scanning telomere units in {fasta} ...")
        units = scan_telomere_units(fasta)
        unit_bed = outdir / f"{prefix}.tel_unit.bed"
        write_unit_bed(units, unit_bed)
        print(f"Wrote {len(units)} deduplicated units -> {unit_bed}")

    print("Computing bin density ...")
    bin_rows = compute_bin_density(units, chrom_lengths, args.bin_size)
    write_density_table(bin_rows, outdir / f"{prefix}_tel.{args.bin_size}.result")
    write_density_detail(
        bin_rows, outdir / f"{prefix}_tel.{args.bin_size}.detail.tsv"
    )

    print("Calling telomere coordinates ...")
    calls = call_all_telomeres(
        bin_rows,
        units,
        chrom_lengths,
        args.density_threshold,
        args.max_unit_gap,
        args.max_dip_bins,
    )
    coords_path = outdir / f"{prefix}.telomere_coords.tsv"
    write_coords(calls, coords_path)
    print(f"Wrote {len(calls)} telomere calls -> {coords_path}")
    for c in calls:
        print(
            f"  {c.chrom}\t{c.end}\t"
            f"tel:{c.tel_start}-{c.tel_end}\t"
            f"len:{c.length_bp}\t{c.dominant_motif}"
        )

    if not args.no_plot:
        try_plot(
            units,
            calls,
            chrom_lengths,
            outdir,
            prefix,
            args.plot_flank,
            args.density_threshold,
            args.bin_size,
        )



# ========== Read coverage ==========

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


def parse_bam_arg(spec: str) -> Tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"BAM spec must be label=path, got: {spec}")
    label, path = spec.split("=", 1)
    label = label.strip()
    path = Path(path.strip())
    if not label:
        raise ValueError(f"Empty BAM label in: {spec}")
    if not path.is_file():
        raise FileNotFoundError(f"BAM not found: {path}")
    return label, path


def load_telomeres(coords_path: Path) -> List[dict]:
    with coords_path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"No telomere rows in {coords_path}")
    required = {"chromosome", "end", "telomere_start", "telomere_end"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Missing columns in {coords_path}: {sorted(missing)}")
    return rows


def make_bins(telomeres: Sequence[dict], bin_size: int) -> List[dict]:
    bins: List[dict] = []
    for tel in telomeres:
        chrom = tel["chromosome"]
        arm = tel["end"]
        t0 = int(tel["telomere_start"])
        t1 = int(tel["telomere_end"])
        if t1 <= t0:
            raise ValueError(f"Invalid telomere span on {chrom} {arm}: {t0}-{t1}")

        grid_start = (t0 // bin_size) * bin_size
        grid_end = ((t1 + bin_size - 1) // bin_size) * bin_size
        for bs in range(grid_start, grid_end, bin_size):
            be = bs + bin_size
            bins.append(
                {
                    "chromosome": chrom,
                    "end": arm,
                    "bin_start": bs,
                    "bin_end": be,
                    "telomere_start": t0,
                    "telomere_end": t1,
                    "dominant_motif": tel.get("dominant_motif", ""),
                    "length_bp": tel.get("length_bp", ""),
                }
            )
    return bins


def write_bins_bed(bins: Sequence[dict], bed_path: Path) -> None:
    with bed_path.open("w") as fh:
        for i, b in enumerate(bins):
            name = f"{b['chromosome']}_{b['end']}_{b['bin_start']}"
            fh.write(
                f"{b['chromosome']}\t{b['bin_start']}\t{b['bin_end']}\t{name}\n"
            )


def _run_bedtools_coverage_mode(
    bedtools: str,
    bed_path: Path,
    bam_path: Path,
    mode: str,
) -> Dict[Tuple[str, int, int], float]:
    """Run bedtools coverage with a single output mode (-counts or -mean)."""
    if mode not in {"counts", "mean"}:
        raise ValueError(f"Unsupported bedtools coverage mode: {mode}")

    cmd = [
        bedtools,
        "coverage",
        "-a",
        str(bed_path),
        "-b",
        str(bam_path),
        f"-{mode}",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bedtools coverage -{mode} failed for {bam_path}:\n{proc.stderr.strip()}"
        )

    values: Dict[Tuple[str, int, int], float] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise RuntimeError(f"Unexpected bedtools output: {line}")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        values[(chrom, start, end)] = float(fields[-1])
    return values


def run_bedtools_coverage(
    bedtools: str, bed_path: Path, bam_path: Path
) -> Dict[Tuple[str, int, int], Tuple[int, float]]:
    if shutil.which(bedtools) is None:
        raise RuntimeError(f"bedtools not found: {bedtools}")

    # bedtools v2.31+ treats -counts and -mean as mutually exclusive.
    counts = _run_bedtools_coverage_mode(bedtools, bed_path, bam_path, "counts")
    means = _run_bedtools_coverage_mode(bedtools, bed_path, bam_path, "mean")

    coverage: Dict[Tuple[str, int, int], Tuple[int, float]] = {}
    for key in counts:
        coverage[key] = (int(counts[key]), means.get(key, 0.0))
    for key in means:
        if key not in coverage:
            coverage[key] = (0, means[key])
    return coverage


def write_coverage_tsv(
    out_path: Path,
    bins: Sequence[dict],
    bam_labels: Sequence[str],
    coverage_by_bam: Dict[str, Dict[Tuple[str, int, int], Tuple[int, float]]],
) -> None:
    base_cols = [
        "chromosome",
        "end",
        "bin_start",
        "bin_end",
        "telomere_start",
        "telomere_end",
        "dominant_motif",
        "length_bp",
    ]
    bam_cols: List[str] = []
    for label in bam_labels:
        bam_cols.extend([f"{label}_read_count", f"{label}_mean_depth"])

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(base_cols + bam_cols)
        for b in bins:
            key = (b["chromosome"], b["bin_start"], b["bin_end"])
            row = [
                b["chromosome"],
                b["end"],
                b["bin_start"],
                b["bin_end"],
                b["telomere_start"],
                b["telomere_end"],
                b["dominant_motif"],
                b["length_bp"],
            ]
            for label in bam_labels:
                count, depth = coverage_by_bam[label].get(key, (0, 0.0))
                row.extend([count, f"{depth:.4f}"])
            writer.writerow(row)


def write_summary_tsv(
    out_path: Path,
    telomeres: Sequence[dict],
    bins: Sequence[dict],
    bam_labels: Sequence[str],
    coverage_by_bam: Dict[str, Dict[Tuple[str, int, int], Tuple[int, float]]],
) -> None:
    cols = [
        "chromosome",
        "end",
        "telomere_start",
        "telomere_end",
        "length_bp",
        "dominant_motif",
        "n_bins",
    ]
    for label in bam_labels:
        cols.extend(
            [
                f"{label}_total_reads",
                f"{label}_mean_depth",
                f"{label}_covered_bins",
                f"{label}_covered_fraction",
            ]
        )

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(cols)
        for tel in telomeres:
            chrom = tel["chromosome"]
            arm = tel["end"]
            t0 = int(tel["telomere_start"])
            t1 = int(tel["telomere_end"])
            tel_bins = [
                b
                for b in bins
                if b["chromosome"] == chrom
                and b["end"] == arm
                and b["telomere_start"] == t0
                and b["telomere_end"] == t1
            ]
            row = [
                chrom,
                arm,
                t0,
                t1,
                tel.get("length_bp", t1 - t0),
                tel.get("dominant_motif", ""),
                len(tel_bins),
            ]
            for label in bam_labels:
                cov = coverage_by_bam[label]
                total_reads = 0
                depth_sum = 0.0
                covered = 0
                for b in tel_bins:
                    key = (b["chromosome"], b["bin_start"], b["bin_end"])
                    count, depth = cov.get(key, (0, 0.0))
                    total_reads += count
                    depth_sum += depth
                    if count > 0:
                        covered += 1
                n = len(tel_bins) or 1
                row.extend(
                    [
                        total_reads,
                        f"{depth_sum / n:.4f}",
                        covered,
                        f"{covered / n:.4f}",
                    ]
                )
            writer.writerow(row)


def main_coverage(args: argparse.Namespace) -> int:
    coords_path = Path(args.coords)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not coords_path.is_file():
        print(f"ERROR: coords not found: {coords_path}", file=sys.stderr)
        return 1

    bam_inputs = [parse_bam_arg(spec) for spec in args.bam]
    telomeres = load_telomeres(coords_path)
    bins = make_bins(telomeres, args.bin_size)

    bed_path = outdir / f"{args.prefix}.telomere_bins.{args.bin_size}bp.bed"
    write_bins_bed(bins, bed_path)
    print(f"Wrote {len(bins)} bins -> {bed_path}")

    coverage_by_bam: Dict[str, Dict[Tuple[str, int, int], Tuple[int, float]]] = {}
    for label, bam_path in bam_inputs:
        bai = Path(f"{bam_path}.bai")
        if not bai.is_file():
            print(
                f"ERROR: BAM index missing: {bai} (run samtools index first)",
                file=sys.stderr,
            )
            return 1
        print(f"Coverage [{label}] <- {bam_path}")
        coverage_by_bam[label] = run_bedtools_coverage(
            args.bedtools, bed_path, bam_path
        )

    labels = [label for label, _ in bam_inputs]
    detail_path = outdir / f"{args.prefix}.telomere_coverage.{args.bin_size}bp.tsv"
    summary_path = outdir / f"{args.prefix}.telomere_coverage.{args.bin_size}bp.summary.tsv"

    write_coverage_tsv(detail_path, bins, labels, coverage_by_bam)
    write_summary_tsv(summary_path, telomeres, bins, labels, coverage_by_bam)

    print(f"Wrote per-bin coverage -> {detail_path}")
    print(f"Wrote telomere summary -> {summary_path}")
    return 0




# ========== Read alignments ==========

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple



@dataclass
class AlignmentRecord:
    chromosome: str
    end: str
    platform: str
    read_name: str
    flag: int
    mapq: int
    strand: str
    ref_start: int
    ref_end: int
    spans_subtelomere: bool
    blocks: List[Tuple[int, int]] = field(default_factory=list)

    def overlap_bp(self, win_start: int, win_end: int) -> int:
        o0 = max(self.ref_start, win_start)
        o1 = min(self.ref_end, win_end)
        return max(0, o1 - o0)


def parse_bam_arg(spec: str) -> Tuple[str, Path]:
    if "=" not in spec:
        raise ValueError(f"BAM spec must be label=path, got: {spec}")
    label, path = spec.split("=", 1)
    label = label.strip()
    path = Path(path.strip())
    if not label:
        raise ValueError(f"Empty BAM label in: {spec}")
    if not path.is_file():
        raise FileNotFoundError(f"BAM not found: {path}")
    return label, path


def load_telomere_rows(coords_path: Path) -> List[dict]:
    with coords_path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in {coords_path}")
    required = {"chromosome", "end", "telomere_start", "telomere_end"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Missing columns in {coords_path}: {sorted(missing)}")
    return rows


def fetch_window(
    chrom: str, end: str, chrom_len: int, flank_bp: int
) -> Tuple[int, int]:
    """0-based half-open [start, end) interval for pysam.fetch."""
    flank_bp = min(flank_bp, chrom_len)
    if end == "5p":
        return 0, flank_bp
    return max(0, chrom_len - flank_bp), chrom_len


def spans_subtelomere_to_telomere(
    ref_start: int,
    ref_end: int,
    tel_start: int,
    tel_end: int,
    end: str,
) -> bool:
    """
    True when the alignment overlaps the telomere call and extends into
    subtelomere (centromere-proximal side of the telomere array).
    """
    if ref_end <= tel_start or ref_start >= tel_end:
        return False
    if end == "5p":
        return ref_end > tel_end
    if end == "3p":
        return ref_start < tel_start
    return False


def encode_blocks(blocks: Sequence[Tuple[int, int]]) -> str:
    return ";".join(f"{a}-{b}" for a, b in blocks)


def decode_blocks(spec: str) -> List[Tuple[int, int]]:
    if not spec:
        return []
    out: List[Tuple[int, int]] = []
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        a, b = part.split("-", 1)
        out.append((int(a), int(b)))
    return out


def alignment_from_pysam(
    read,
    chrom: str,
    end: str,
    platform: str,
    tel_start: int,
    tel_end: int,
    win_start: int,
    win_end: int,
    min_overlap: int,
) -> Optional[AlignmentRecord]:
    if read.is_unmapped:
        return None
    if read.reference_name != chrom:
        return None

    ref_start = int(read.reference_start)
    ref_end = int(read.reference_end)
    overlap = min(ref_end, win_end) - max(ref_start, win_start)
    if overlap < min_overlap:
        return None

    blocks = [(int(a), int(b)) for a, b in read.get_blocks()]
    if not blocks:
        return None

    strand = "-" if read.is_reverse else "+"
    spans = spans_subtelomere_to_telomere(
        ref_start, ref_end, tel_start, tel_end, end
    )
    return AlignmentRecord(
        chromosome=chrom,
        end=end,
        platform=platform,
        read_name=read.query_name,
        flag=int(read.flag),
        mapq=int(read.mapping_quality) if read.mapping_quality is not None else 0,
        strand=strand,
        ref_start=ref_start,
        ref_end=ref_end,
        spans_subtelomere=spans,
        blocks=blocks,
    )


def dedupe_best_per_read(records: Sequence[AlignmentRecord]) -> List[AlignmentRecord]:
    best: Dict[str, AlignmentRecord] = {}
    for rec in records:
        key = rec.read_name
        prev = best.get(key)
        if prev is None:
            best[key] = rec
            continue
        prev_score = (prev.spans_subtelomere, prev.ref_end - prev.ref_start)
        score = (rec.spans_subtelomere, rec.ref_end - rec.ref_start)
        if score > prev_score:
            best[key] = rec
    return list(best.values())


def sort_for_export(records: Sequence[AlignmentRecord]) -> List[AlignmentRecord]:
    return sorted(
        records,
        key=lambda r: (
            0 if r.spans_subtelomere else 1,
            -(r.ref_end - r.ref_start),
            r.ref_start,
            r.read_name,
        ),
    )


def cap_reads(
    records: Sequence[AlignmentRecord], max_reads: int
) -> List[AlignmentRecord]:
    if len(records) <= max_reads:
        return list(records)
    spanning = [r for r in records if r.spans_subtelomere]
    other = [r for r in records if not r.spans_subtelomere]
    out = spanning[:max_reads]
    if len(out) < max_reads:
        out.extend(other[: max_reads - len(out)])
    return out


def extract_arm(
    bam_path: Path,
    platform: str,
    tel_row: dict,
    chrom_len: int,
    flank_bp: int,
    min_mapq: int,
    min_overlap: int,
    max_reads: int,
    include_secondary: bool,
    include_supplementary: bool,
) -> List[AlignmentRecord]:
    import pysam

    chrom = tel_row["chromosome"]
    end = tel_row["end"]
    tel_start = int(tel_row["telomere_start"])
    tel_end = int(tel_row["telomere_end"])
    win_start, win_end = fetch_window(chrom, end, chrom_len, flank_bp)

    raw: List[AlignmentRecord] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if chrom not in bam.references:
            raise ValueError(f"{chrom} not in BAM {bam_path}")
        for read in bam.fetch(chrom, win_start, win_end):
            mq = 0 if read.mapping_quality is None else int(read.mapping_quality)
            if mq < min_mapq:
                continue
            if read.is_secondary and not include_secondary:
                continue
            if read.is_supplementary and not include_supplementary:
                continue
            rec = alignment_from_pysam(
                read,
                chrom,
                end,
                platform,
                tel_start,
                tel_end,
                win_start,
                win_end,
                min_overlap,
            )
            if rec is not None:
                raw.append(rec)

    deduped = dedupe_best_per_read(raw)
    capped = cap_reads(sort_for_export(deduped), max_reads)
    return capped


def write_alignments_tsv(records: Sequence[AlignmentRecord], out_path: Path) -> None:
    fields = [
        "chromosome",
        "end",
        "platform",
        "read_name",
        "flag",
        "mapq",
        "strand",
        "ref_start",
        "ref_end",
        "spans_subtelomere",
        "blocks",
    ]
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rec in records:
            writer.writerow(
                {
                    "chromosome": rec.chromosome,
                    "end": rec.end,
                    "platform": rec.platform,
                    "read_name": rec.read_name,
                    "flag": rec.flag,
                    "mapq": rec.mapq,
                    "strand": rec.strand,
                    "ref_start": rec.ref_start,
                    "ref_end": rec.ref_end,
                    "spans_subtelomere": int(rec.spans_subtelomere),
                    "blocks": encode_blocks(rec.blocks),
                }
            )


def load_alignments_tsv(path: Path) -> List[AlignmentRecord]:
    csv.field_size_limit(sys.maxsize)
    records: List[AlignmentRecord] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            records.append(
                AlignmentRecord(
                    chromosome=row["chromosome"],
                    end=row["end"],
                    platform=row["platform"],
                    read_name=row["read_name"],
                    flag=int(row["flag"]),
                    mapq=int(row["mapq"]),
                    strand=row["strand"],
                    ref_start=int(row["ref_start"]),
                    ref_end=int(row["ref_end"]),
                    spans_subtelomere=bool(int(row["spans_subtelomere"])),
                    blocks=decode_blocks(row["blocks"]),
                )
            )
    return records


def index_by_arm_platform(
    records: Sequence[AlignmentRecord],
) -> Dict[Tuple[str, str, str], List[AlignmentRecord]]:
    out: Dict[Tuple[str, str, str], List[AlignmentRecord]] = {}
    for rec in records:
        key = (rec.chromosome, rec.end, rec.platform)
        out.setdefault(key, []).append(rec)
    return out


def block_to_plot_x(
    bstart: int, bend: int, end: str, chrom_len: int
) -> Tuple[float, float]:
    """Map 0-based genomic block to plot x (distance from chromosome end)."""
    if end == "5p":
        return float(bstart), float(bend)
    return float(chrom_len - bend), float(chrom_len - bstart)


def assign_plot_lanes(
    records: Sequence[AlignmentRecord],
    end: str,
    chrom_len: int,
    flank_bp: int,
) -> Dict[str, int]:
    """Greedy lane packing in plot-x coordinates (IGV-style)."""

    def plot_extent(rec: AlignmentRecord) -> Tuple[float, float]:
        xs: List[float] = []
        xe: List[float] = []
        for b0, b1 in rec.blocks:
            xl, xr = block_to_plot_x(b0, b1, end, chrom_len)
            xs.append(xl)
            xe.append(xr)
        return min(xs), max(xe)

    ordered = sorted(
        records,
        key=lambda r: (
            0 if r.spans_subtelomere else 1,
            plot_extent(r)[0],
            r.read_name,
        ),
    )
    lane_ends: List[float] = []
    lanes: Dict[str, int] = {}
    gap = 5.0
    for rec in ordered:
        x0, x1 = plot_extent(rec)
        x0 = max(0.0, x0)
        x1 = min(float(flank_bp), x1)
        if x1 <= x0:
            continue
        placed = False
        for i, end_x in enumerate(lane_ends):
            if x0 >= end_x + gap:
                lane_ends[i] = x1
                lanes[rec.read_name] = i
                placed = True
                break
        if not placed:
            lanes[rec.read_name] = len(lane_ends)
            lane_ends.append(x1)
    return lanes


def main_alignments(args: argparse.Namespace) -> int:
    coords_path = Path(args.coords)
    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        import pysam  # noqa: F401
    except ImportError:
        print(
            "ERROR: pysam is required. Example: conda install -c bioconda pysam",
            file=sys.stderr,
        )
        return 1

    for path, label in ((coords_path, "coords"), (fasta_path, "FASTA")):
        if not path.is_file():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 1

    bam_specs = [parse_bam_arg(s) for s in args.bam]
    tel_rows = load_telomere_rows(coords_path)
    fai = resolve_fai(fasta_path, args.fai)
    chrom_lengths = load_chrom_lengths(fai)

    all_records: List[AlignmentRecord] = []
    for platform, bam_path in bam_specs:
        print(f"Scanning {platform}: {bam_path}")
        for tel_row in tel_rows:
            chrom = tel_row["chromosome"]
            arm = tel_row["end"]
            chrom_len = chrom_lengths[chrom]
            arm_recs = extract_arm(
                bam_path,
                platform,
                tel_row,
                chrom_len,
                args.plot_flank,
                args.min_mapq,
                args.min_overlap,
                args.max_reads_per_arm,
                args.include_secondary,
                args.include_supplementary,
            )
            n_span = sum(1 for r in arm_recs if r.spans_subtelomere)
            print(
                f"  {chrom} {arm}: {len(arm_recs)} alignments "
                f"({n_span} span subtelomere→telomere)"
            )
            all_records.extend(arm_recs)

    out_tsv = outdir / f"{args.prefix}.telomere_alignments.tsv"
    write_alignments_tsv(all_records, out_tsv)
    print(f"Wrote {len(all_records)} rows -> {out_tsv}")
    return 0




# ========== Plotting ==========

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Distinct from MOTIF_COLORS (CCCTAAA blue / TTTAGGG red) used in the panel above.
# Vertical layout: top span track, motif density, read-depth bars (shorter panels).

def _configure_plot_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.grid": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _detect_platform_labels(header: Sequence[str]) -> List[str]:
    labels: List[str] = []
    for col in header:
        if col.endswith("_mean_depth"):
            labels.append(col[: -len("_mean_depth")])
    if not labels:
        raise ValueError("No *_mean_depth columns found in coverage TSV")
    return labels


def load_coverage(path: Path) -> Tuple[List[str], List[dict]]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid TSV: {path}")
        platforms = _detect_platform_labels(reader.fieldnames)
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return platforms, rows


def load_summary(path: Path) -> Tuple[List[str], List[dict]]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty or invalid summary TSV: {path}")
        platforms = _detect_platform_labels(reader.fieldnames)
        rows = list(reader)
    return platforms, rows


def load_calls_from_coords(path: Path) -> List[TelomereCall]:
    calls: List[TelomereCall] = []
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            calls.append(
                TelomereCall(
                    chrom=row["chromosome"],
                    end=row["end"],
                    bin_start=int(row["bin_region_start"]),
                    bin_end=int(row["bin_region_end"]),
                    tel_start=int(row["telomere_start"]),
                    tel_end=int(row["telomere_end"]),
                    length_bp=int(row["length_bp"]),
                    dominant_motif=row["dominant_motif"],
                    n_units=int(row["n_units"]),
                    mean_bin_density=float(row["mean_bin_density"]),
                    max_bin_density=float(row["max_bin_density"]),
                )
            )
    return calls


def _bin_center(row: dict) -> float:
    return (int(row["bin_start"]) + int(row["bin_end"])) / 2.0


def _x_from_end(row: dict, end: str, chrom_len: int) -> float:
    center = _bin_center(row)
    if end == "5p":
        return center
    return chrom_len - center


def _index_coverage(
    rows: Sequence[dict],
) -> Dict[Tuple[str, str], List[dict]]:
    by_arm: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in rows:
        by_arm[(row["chromosome"], row["end"])].append(row)
    for arm_rows in by_arm.values():
        arm_rows.sort(key=lambda r: int(r["bin_start"]))
    return by_arm


def _draw_combined_end_column(
    fig,
    col_gs,
    chrom: str,
    chrom_len: int,
    end: str,
    call: Optional[TelomereCall],
    panel_bins: Sequence[Tuple[float, float, float]],
    cov_rows: Sequence[dict],
    platforms: Sequence[str],
    flank_bp: int,
    bin_size: int,
    density_threshold: float,
) -> None:
    inner = col_gs.subgridspec(3, 1, height_ratios=ENDS_HEIGHT_RATIOS, hspace=0.06)
    ax_top = fig.add_subplot(inner[0])
    ax_mid = fig.add_subplot(inner[1], sharex=ax_top)
    ax_cov = fig.add_subplot(inner[2], sharex=ax_top)

    ax_top.set_xlim(0, flank_bp)
    if end == "3p":
        ax_top.invert_xaxis()

    ax_top.set_ylim(0, 1)
    ax_top.set_yticks([])
    span = _tel_span_plot_coords(call, chrom_len, end, flank_bp)
    if span is not None and call is not None:
        x0, x1 = span
        ax_top.plot(
            [x0, x1],
            [0.38, 0.38],
            color="#C73E1D",
            solid_capstyle="butt",
            linewidth=3.5,
            clip_on=False,
            zorder=3,
        )
        ax_top.text(
            (x0 + x1) / 2,
            0.58,
            f"{_format_bp(call.tel_start)} – {_format_bp(call.tel_end)} bp",
            ha="center",
            va="bottom",
            fontsize=7,
            fontfamily="Liberation Sans",
            zorder=4,
        )
    elif call is None:
        ax_top.text(
            flank_bp / 2,
            0.45,
            "No telomere call",
            ha="center",
            va="center",
            fontsize=7,
            fontfamily="Liberation Sans",
            color="#888888",
        )

    ax_top.tick_params(
        axis="x",
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        length=4,
        width=0.8,
        labelsize=7,
        direction="in",
    )
    for lbl in ax_top.get_xticklabels():
        lbl.set_fontfamily("Liberation Sans")
    ax_top.spines["top"].set_visible(True)
    ax_top.spines["top"].set_linewidth(0.6)
    ax_top.spines["bottom"].set_visible(False)
    ax_top.spines["left"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.set_title(f"{chrom} {end}", fontsize=8, fontfamily="Liberation Sans", pad=3)

    if not panel_bins:
        ax_mid.text(
            0.5,
            0.5,
            "No bins in window",
            transform=ax_mid.transAxes,
            ha="center",
            fontfamily="Liberation Sans",
        )
    else:
        xs = [b[0] for b in panel_bins]
        ccct = [b[1] for b in panel_bins]
        ttta = [b[2] for b in panel_bins]
        width = bin_size * 0.92

        ax_mid.bar(
            xs,
            ccct,
            width=width,
            label="CCCTAAA",
            color=MOTIF_COLORS["CCCTAAA"],
            edgecolor="none",
            align="center",
        )
        ax_mid.bar(
            xs,
            ttta,
            width=width,
            bottom=ccct,
            label="TTTAGGG",
            color=MOTIF_COLORS["TTTAGGG"],
            edgecolor="none",
            align="center",
        )
        ax_mid.axhline(
            density_threshold, color="#888888", ls="--", lw=0.9, dashes=(4, 3)
        )

    ax_mid.set_ylim(0, 1.05)
    ax_mid.tick_params(axis="x", labelbottom=False, bottom=False, top=False)
    ax_mid.set_ylabel("Unit density", fontsize=7, fontfamily="Liberation Sans")
    ax_mid.tick_params(axis="y", labelsize=7)
    for lbl in ax_mid.get_yticklabels():
        lbl.set_fontfamily("Liberation Sans")
    legend_loc = "upper left" if end == "3p" else "upper right"
    ax_mid.legend(
        loc=legend_loc,
        frameon=False,
        fontsize=7,
        prop={"family": "Liberation Sans"},
    )
    ax_mid.grid(False)
    for spine in ("top", "right"):
        ax_mid.spines[spine].set_visible(False)

    if cov_rows:
        xs = [_x_from_end(r, end, chrom_len) for r in cov_rows]
        n_plat = len(platforms)
        group_width = bin_size * 0.92
        bar_width = group_width / n_plat
        for i, platform in enumerate(platforms):
            depths = [float(r[f"{platform}_mean_depth"]) for r in cov_rows]
            offset = (i - (n_plat - 1) / 2) * bar_width
            color = PLATFORM_COLORS.get(platform.lower())
            ax_cov.bar(
                [x + offset for x in xs],
                depths,
                width=bar_width * 0.95,
                label=platform.upper(),
                color=color,
                edgecolor="none",
                align="center",
            )
    else:
        ax_cov.text(
            0.5,
            0.5,
            "No coverage bins",
            transform=ax_cov.transAxes,
            ha="center",
            fontfamily="Liberation Sans",
            color="#888888",
        )

    ax_cov.set_xlabel(
        "Distance from telomere end (bp)", fontsize=7, fontfamily="Liberation Sans", labelpad=3
    )
    ax_cov.set_ylabel("Mean read depth", fontsize=7, fontfamily="Liberation Sans")
    ax_cov.tick_params(axis="x", labelbottom=False, bottom=False, top=False)
    ax_cov.tick_params(axis="y", labelsize=7)
    for lbl in ax_cov.get_yticklabels():
        lbl.set_fontfamily("Liberation Sans")
    cov_legend_loc = "upper left" if end == "3p" else "upper right"
    ax_cov.legend(
        loc=cov_legend_loc,
        frameon=False,
        fontsize=7,
        prop={"family": "Liberation Sans"},
    )
    ax_cov.grid(False)
    for spine in ("top", "right"):
        ax_cov.spines[spine].set_visible(False)

    return ax_top, ax_mid, ax_cov


def plot_combined_ends(
    units: Sequence[TelUnit],
    calls: Sequence[TelomereCall],
    chrom_lengths: Dict[str, int],
    coverage_by_arm: Dict[Tuple[str, str], List[dict]],
    platforms: Sequence[str],
    outdir: Path,
    flank_bp: int,
    bin_size: int,
    density_threshold: float,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    calls_by = {(c.chrom, c.end): c for c in calls}

    chroms = sorted(chrom_lengths, key=_natural_key)
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    outer_gs = GridSpec(3, 2, figure=fig, hspace=0.25, wspace=0.35,
                        left=0.08, right=0.96, top=0.97, bottom=0.04)

    for idx, chrom in enumerate(chroms):
        chrom_len = chrom_lengths[chrom]
        row, col = idx // 2, idx % 2
        inner_gs = outer_gs[row, col].subgridspec(1, 2, wspace=0.18)

        for end_col, end in enumerate(("5p", "3p")):
            call = calls_by.get((chrom, end))
            panel_bins = _plot_panel_bins(
                units, chrom, chrom_len, end, flank_bp, bin_size
            )
            cov_rows = coverage_by_arm.get((chrom, end), [])
            _draw_combined_end_column(
                fig,
                inner_gs[end_col],
                chrom,
                chrom_len,
                end,
                call,
                panel_bins,
                cov_rows,
                platforms,
                flank_bp,
                bin_size,
                density_threshold,
            )

    fig.suptitle(
        f"Telomere ends (0–{_format_bp(flank_bp)} bp from each end)",
        fontsize=10,
        fontfamily="Liberation Sans",
        y=0.995,
    )
    fig.savefig(outdir / "all_chromosomes_ends.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_summary_bars(
    summary_rows: Sequence[dict],
    platforms: Sequence[str],
    outdir: Path,
    prefix: str,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [f"{r['chromosome']}_{r['end']}" for r in summary_rows]
    x = np.arange(len(labels))
    n_plat = len(platforms)
    width = 0.8 / n_plat

    fig, ax = plt.subplots(figsize=SUMMARY_FIGSIZE)
    for i, platform in enumerate(platforms):
        depths = [float(r[f"{platform}_mean_depth"]) for r in summary_rows]
        offset = (i - (n_plat - 1) / 2) * width
        color = PLATFORM_COLORS.get(platform.lower())
        ax.bar(
            x + offset,
            depths,
            width=width,
            label=platform.upper(),
            color=color,
            edgecolor="none",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontfamily="Liberation Sans")
    ax.set_ylabel("Mean read depth (100-bp bins)", fontfamily="Liberation Sans")
    ax.set_title("Telomere-region mean read depth per end", fontfamily="Liberation Sans")
    ax.legend(frameon=False, prop={"family": "Liberation Sans"})
    ax.grid(False)
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily("Liberation Sans")
    fig.tight_layout()
    fig.savefig(
        outdir / f"{prefix}_telomere_coverage_summary.pdf", dpi=150, bbox_inches="tight"
    )
    plt.close(fig)


def infer_prefix(input_path: Path) -> str:
    name = input_path.name
    for suffix in (
        ".telomere_coverage.100bp.tsv",
        ".telomere_coverage.tsv",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return input_path.stem


def infer_summary_path(input_path: Path) -> Path:
    return input_path.with_name(f"{input_path.stem}.summary.tsv")


def run_coverage_plots(
    coverage_path: Path,
    coords_path: Path,
    units_path: Path,
    fasta_path: Path,
    outdir: Path,
    *,
    fai: Optional[str] = None,
    prefix: Optional[str] = None,
    summary_path: Optional[Path] = None,
    bin_size: int = 100,
    plot_flank: int = DEFAULT_PLOT_FLANK,
    density_threshold: float = DEFAULT_DENSITY,
) -> int:
    """Render coverage ends figures and optional summary bar chart."""
    for path, label in (
        (coverage_path, "coverage TSV"),
        (coords_path, "coords TSV"),
        (units_path, "tel unit BED"),
        (fasta_path, "FASTA"),
    ):
        if not path.is_file():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 1

    try:
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        print("ERROR: matplotlib is required for plotting.", file=sys.stderr)
        return 1

    _configure_plot_style()
    outdir.mkdir(parents=True, exist_ok=True)
    file_prefix = prefix or infer_prefix(coverage_path)

    fai_path = resolve_fai(fasta_path, fai)
    chrom_lengths = load_chrom_lengths(fai_path)
    units = load_units_bed(units_path)
    calls = load_calls_from_coords(coords_path)
    platforms, cov_rows = load_coverage(coverage_path)
    coverage_by_arm = _index_coverage(cov_rows)

    print(f"Loaded {len(cov_rows)} coverage bins, platforms: {', '.join(platforms)}")
    print(f"Loaded {len(units)} telomere units, {len(calls)} telomere calls")

    plot_combined_ends(
        units,
        calls,
        chrom_lengths,
        coverage_by_arm,
        platforms,
        outdir,
        plot_flank,
        bin_size,
        density_threshold,
    )
    print(f"Wrote combined ends+coverage PDF -> {outdir}/all_chromosomes_ends.pdf")

    summary = summary_path if summary_path else infer_summary_path(coverage_path)
    if summary.is_file():
        _, summary_rows = load_summary(summary)
        plot_summary_bars(summary_rows, platforms, outdir, file_prefix)
        print(
            f"Wrote summary bar chart -> {outdir}/{file_prefix}_telomere_coverage_summary.pdf"
        )
    else:
        print(
            f"Summary TSV not found (skipping bar chart): {summary}",
            file=sys.stderr,
        )

    return 0


import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# outer column: [annotation+coverage | read tracks]

def _platform_read_ylabel(platform: str) -> str:
    label = platform.lower()
    if label == "hifi":
        return "HiFi reads"
    if label == "ont":
        return "ONT reads"
    return f"{platform.upper()} reads"


def _platform_track_color(platform: str) -> str:
    return PLATFORM_COLORS.get(platform.lower(), "#888888")


def _draw_read_track(
    ax,
    records: Sequence[AlignmentRecord],
    end: str,
    chrom_len: int,
    flank_bp: int,
    platform: str,
    lane_height: float,
) -> None:
    ylabel = _platform_read_ylabel(platform)
    records = [r for r in records if r.spans_subtelomere]
    block_color = _platform_track_color(platform)

    if not records:
        ax.text(
            0.5,
            0.5,
            f"No {ylabel}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7,
            fontfamily="Liberation Sans",
            color="#888888",
        )
        ax.set_ylabel(ylabel, fontsize=7, fontfamily="Liberation Sans")
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        return

    lanes = assign_plot_lanes(records, end, chrom_len, flank_bp)
    n_lanes = max(lanes.values()) + 1 if lanes else 1
    ax.set_ylim(-0.5, n_lanes * lane_height + 0.3)
    ax.set_yticks([])

    for rec in records:
        lane = lanes.get(rec.read_name, 0)
        y = lane * lane_height + lane_height * 0.5

        prev_xr: Optional[float] = None
        for b0, b1 in rec.blocks:
            xl, xr = block_to_plot_x(b0, b1, end, chrom_len)
            xl = max(0.0, xl)
            xr = min(float(flank_bp), xr)
            if xr <= xl:
                continue
            ax.plot(
                [xl, xr],
                [y, y],
                color=block_color,
                linewidth=2.0,
                alpha=0.9,
                solid_capstyle="butt",
                zorder=3,
            )
            if prev_xr is not None and xl > prev_xr + 1:
                ax.plot(
                    [prev_xr, xl],
                    [y, y],
                    color=block_color,
                    linewidth=0.6,
                    alpha=0.4,
                    zorder=1,
                )
            prev_xr = xr

    ax.set_ylabel(ylabel, fontsize=7, fontfamily="Liberation Sans")
    ax.tick_params(axis="x", labelbottom=False, bottom=False, top=False)
    ax.tick_params(axis="y", labelsize=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(False)


def _draw_end_column_with_reads(
    fig,
    col_gs,
    chrom: str,
    chrom_len: int,
    end: str,
    call: Optional[TelomereCall],
    panel_bins,
    cov_rows,
    platforms: Sequence[str],
    align_by_plat: Dict[str, List[AlignmentRecord]],
    flank_bp: int,
    bin_size: int,
    density_threshold: float,
    lane_height: float,
) -> None:
    outer = col_gs.subgridspec(
        2, 1, height_ratios=COLUMN_OUTER_RATIOS, hspace=0.08
    )
    ax_top, _, ax_cov = _draw_combined_end_column(
        fig,
        outer[0],
        chrom,
        chrom_len,
        end,
        call,
        panel_bins,
        cov_rows,
        platforms,
        flank_bp,
        bin_size,
        density_threshold,
    )

    n_plat = len(platforms)
    read_gs = outer[1].subgridspec(
        n_plat,
        1,
        height_ratios=[READ_TRACK_ROW_RATIO] * n_plat,
        hspace=0.18,
    )
    for i, platform in enumerate(platforms):
        ax_reads = fig.add_subplot(read_gs[i], sharex=ax_top)
        recs = align_by_plat.get(platform, [])
        _draw_read_track(
            ax_reads,
            recs,
            end,
            chrom_len,
            flank_bp,
            platform,
            lane_height,
        )
        if i < n_plat - 1:
            ax_reads.tick_params(axis="x", labelbottom=False)

    ax_cov.set_xlabel(
        "Distance from telomere end (bp)",
        fontsize=7,
        fontfamily="Liberation Sans",
        labelpad=3,
    )


def plot_ends_with_reads(
    units: Sequence[TelUnit],
    calls: Sequence[TelomereCall],
    chrom_lengths: Dict[str, int],
    coverage_by_arm: Dict[Tuple[str, str], List[dict]],
    platforms: Sequence[str],
    align_index: Dict[Tuple[str, str, str], List[AlignmentRecord]],
    outdir: Path,
    flank_bp: int,
    bin_size: int,
    density_threshold: float,
    suffix: str,
    lane_height: float,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    calls_by = {(c.chrom, c.end): c for c in calls}
    outdir.mkdir(parents=True, exist_ok=True)

    chroms = sorted(chrom_lengths, key=_natural_key)
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    outer_gs = GridSpec(3, 2, figure=fig, hspace=0.25, wspace=0.35,
                        left=0.08, right=0.96, top=0.97, bottom=0.04)

    for idx, chrom in enumerate(chroms):
        chrom_len = chrom_lengths[chrom]
        row, col = idx // 2, idx % 2
        inner_gs = outer_gs[row, col].subgridspec(1, 2, wspace=0.18)

        for end_col, end in enumerate(("5p", "3p")):
            call = calls_by.get((chrom, end))
            panel_bins = _plot_panel_bins(
                units, chrom, chrom_len, end, flank_bp, bin_size
            )
            cov_rows = coverage_by_arm.get((chrom, end), [])
            align_by_plat = {
                plat: align_index.get((chrom, end, plat), [])
                for plat in platforms
            }
            _draw_end_column_with_reads(
                fig,
                inner_gs[end_col],
                chrom,
                chrom_len,
                end,
                call,
                panel_bins,
                cov_rows,
                platforms,
                align_by_plat,
                flank_bp,
                bin_size,
                density_threshold,
                lane_height,
            )

    fig.suptitle(
        f"Telomere ends + read alignments (0–{_format_bp(flank_bp)} bp) — {suffix.lstrip('_ends_')}",
        fontsize=10,
        fontfamily="Liberation Sans",
        y=0.995,
    )
    out_name = f"all_chromosomes{suffix}.pdf"
    fig.savefig(outdir / out_name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outdir / out_name}")


def run_reads_plots(
    alignments_path: Path,
    coords_path: Path,
    units_path: Path,
    fasta_path: Path,
    outdir: Path,
    *,
    fai: Optional[str] = None,
    coverage_path: Optional[Path] = None,
    bin_size: int = 100,
    plot_flank: int = DEFAULT_PLOT_FLANK,
    density_threshold: float = DEFAULT_DENSITY,
    suffix: str = "_ends_reads",
    lane_height: float = 0.85,
) -> int:
    """Render ends figures with spanning-read alignment tracks."""
    for path, label in (
        (alignments_path, "alignments TSV"),
        (coords_path, "coords TSV"),
        (units_path, "tel unit BED"),
        (fasta_path, "FASTA"),
    ):
        if not path.is_file():
            print(f"ERROR: {label} not found: {path}", file=sys.stderr)
            return 1

    try:
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        print("ERROR: matplotlib is required.", file=sys.stderr)
        return 1

    _configure_plot_style()

    fai_path = resolve_fai(fasta_path, fai)
    chrom_lengths = load_chrom_lengths(fai_path)
    units = load_units_bed(units_path)
    calls = load_calls_from_coords(coords_path)
    align_records = load_alignments_tsv(alignments_path)
    align_index = index_by_arm_platform(align_records)

    platforms = sorted({r.platform for r in align_records})
    if not platforms:
        print("ERROR: no platforms in alignment TSV", file=sys.stderr)
        return 1

    coverage_by_arm: Dict[Tuple[str, str], List[dict]] = {}
    if coverage_path is not None:
        if not coverage_path.is_file():
            print(f"ERROR: coverage TSV not found: {coverage_path}", file=sys.stderr)
            return 1
        cov_platforms, cov_rows = load_coverage(coverage_path)
        for p in cov_platforms:
            if p not in platforms:
                platforms.append(p)
        platforms = sorted(set(platforms))
        coverage_by_arm = _index_coverage(cov_rows)

    outdir.mkdir(parents=True, exist_ok=True)
    plot_ends_with_reads(
        units,
        calls,
        chrom_lengths,
        coverage_by_arm,
        platforms,
        align_index,
        outdir,
        plot_flank,
        bin_size,
        density_threshold,
        suffix,
        lane_height,
    )
    return 0




# ========== CLI ==========


PLATFORM_COLORS = {"ont": "#43AA8B", "hifi": "#F8961E"}
ENDS_HEIGHT_RATIOS = [0.70, 2.1, 1.20]
ENDS_FIGSIZE = (14, 6.8)
SUMMARY_FIGSIZE = (10, 3.6)
COLUMN_OUTER_RATIOS = [5.4, 2.6]
READ_TRACK_ROW_RATIO = 1.0
ENDS_WITH_READS_FIGSIZE = (14, 10.5)


@dataclass
class TelomerePaths:
    workdir: Path
    prefix: str
    bin_size: int = 100
    plot_flank: int = DEFAULT_PLOT_FLANK

    @property
    def coords(self) -> Path:
        return self.workdir / f"{self.prefix}.telomere_coords.tsv"

    @property
    def units(self) -> Path:
        return self.workdir / f"{self.prefix}.tel_unit.bed"

    @property
    def coverage(self) -> Path:
        return self.workdir / f"{self.prefix}.telomere_coverage.{self.bin_size}bp.tsv"

    @property
    def coverage_summary(self) -> Path:
        return self.workdir / f"{self.prefix}.telomere_coverage.{self.bin_size}bp.summary.tsv"

    @property
    def alignments(self) -> Path:
        return self.workdir / f"{self.prefix}.telomere_alignments.tsv"

    @property
    def plot_dir(self) -> Path:
        return self.workdir / "plots"


def cmd_find_tel(args: argparse.Namespace) -> int:
    out = Path(args.workdir)
    out.mkdir(parents=True, exist_ok=True)
    units = out / f"{args.prefix}.tel_unit.bed"
    reuse_path = None
    if units.is_file():
        reuse_path = str(units)
    elif args.reuse_units:
        print(
            f"NOTE: {units} not found; scanning telomere units "
            "(--reuse-units ignored)"
        )

    ns = argparse.Namespace(
        fasta=args.fasta_asm,
        fai=args.fai_asm,
        outdir=str(out),
        prefix=args.prefix,
        reuse_units=reuse_path,
        bin_size=args.bin_size,
        density_threshold=args.density_threshold,
        max_unit_gap=DEFAULT_MAX_UNIT_GAP,
        max_dip_bins=DEFAULT_MAX_DIP_BINS,
        plot_flank=args.plot_flank,
        no_plot=True,
    )
    run_pipeline(ns)
    return 0


def cmd_coverage(args: argparse.Namespace) -> int:
    paths = TelomerePaths(Path(args.workdir), args.prefix, args.bin_size)
    if not paths.coords.is_file():
        print(f"ERROR: run find-tel first: {paths.coords}", file=sys.stderr)
        return 1
    ns = argparse.Namespace(
        coords=str(paths.coords),
        bam=[f"ont={args.ont_bam}", f"hifi={args.hifi_bam}"],
        outdir=str(paths.workdir),
        prefix=args.prefix,
        bin_size=args.bin_size,
        bedtools=args.bedtools,
    )
    return main_coverage(ns)


def cmd_alignments(args: argparse.Namespace) -> int:
    paths = TelomerePaths(Path(args.workdir), args.prefix, args.bin_size, args.plot_flank)
    if not paths.coords.is_file():
        print(f"ERROR: run find-tel first: {paths.coords}", file=sys.stderr)
        return 1
    try:
        import pysam  # noqa: F401
    except ImportError:
        print("ERROR: pysam required. conda install -c bioconda pysam", file=sys.stderr)
        return 1
    ns = argparse.Namespace(
        coords=str(paths.coords),
        fasta=args.fasta_plot,
        fai=args.fai_plot,
        bam=[f"ont={args.ont_bam}", f"hifi={args.hifi_bam}"],
        outdir=str(paths.workdir),
        prefix=args.prefix,
        plot_flank=args.plot_flank,
        min_mapq=args.min_mapq,
        min_overlap=args.min_overlap,
        max_reads_per_arm=args.max_reads_per_arm,
        include_secondary=False,
        include_supplementary=False,
    )
    return main_alignments(ns)


def cmd_plot(args: argparse.Namespace) -> int:
    paths = TelomerePaths(Path(args.workdir), args.prefix, args.bin_size, args.plot_flank)
    paths.plot_dir.mkdir(parents=True, exist_ok=True)
    if not paths.coords.is_file() or not paths.units.is_file():
        print("ERROR: coords/units missing; run find-tel first.", file=sys.stderr)
        return 1
    try:
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        print("ERROR: matplotlib required.", file=sys.stderr)
        return 1

    status = 0
    calls = load_calls_from_coords(paths.coords)
    print("--- telomere tract lengths ---")
    plot_telomere_lengths(calls, paths.plot_dir, args.prefix)

    if paths.coverage.is_file():
        print("--- ends + coverage ---")
        status |= run_coverage_plots(
            paths.coverage,
            paths.coords,
            paths.units,
            Path(args.fasta_plot),
            paths.plot_dir,
            fai=args.fai_plot,
            prefix=args.prefix,
            summary_path=paths.coverage_summary,
            bin_size=args.bin_size,
            plot_flank=args.plot_flank,
            density_threshold=args.density_threshold,
        )
    else:
        print(f"Skip coverage plots: {paths.coverage}", file=sys.stderr)

    if paths.alignments.is_file():
        print("--- ends + spanning reads ---")
        status |= run_reads_plots(
            paths.alignments,
            paths.coords,
            paths.units,
            Path(args.fasta_plot),
            paths.plot_dir,
            fai=args.fai_plot,
            coverage_path=paths.coverage if paths.coverage.is_file() else None,
            bin_size=args.bin_size,
            plot_flank=args.plot_flank,
            density_threshold=args.density_threshold,
        )
    else:
        print(f"Skip read plots: {paths.alignments}", file=sys.stderr)

    print(f"Figures -> {paths.plot_dir}/")
    return status


def cmd_all(args: argparse.Namespace) -> int:
    for fn in (cmd_find_tel, cmd_coverage, cmd_alignments, cmd_plot):
        rc = fn(args)
        if rc:
            return rc
    return 0


def _add_shared_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--workdir", type=Path, default=Path("."))
    p.add_argument("--prefix", default="GifuT2T")
    p.add_argument("--bin-size", type=int, default=100)
    p.add_argument("--plot-flank", type=int, default=DEFAULT_PLOT_FLANK)
    p.add_argument("--density-threshold", type=float, default=DEFAULT_DENSITY)
    p.add_argument(
        "--fasta-asm",
        default="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta",
    )
    p.add_argument("--fai-asm", default=None)
    p.add_argument(
        "--fasta-plot",
        default="{PROJ_04_LOTUS_GENOME}/00_assembly_result/GifuT2T/Lotus_GifuT2T_v1.0.fasta",
    )
    p.add_argument("--fai-plot", default=None)
    p.add_argument(
        "--ont-bam",
        default="{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/08.quality/01.mapping/01.ont/ont_q30.bam",
    )
    p.add_argument(
        "--hifi-bam",
        default="{PROJ_04_LOTUS_GENOME}/02_assembly/Gifu/08.quality/01.mapping/02.hifi/ccs_q30.bam",
    )
    p.add_argument(
        "--bedtools",
        default="{CONDA_PY3}/bin/bedtools",
    )
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--min-overlap", type=int, default=50)
    p.add_argument("--max-reads-per-arm", type=int, default=2500)
    p.add_argument("--reuse-units", action="store_true")


def build_main_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Gifu telomere analysis")
    sub = p.add_subparsers(dest="command", required=True)
    steps = (
        ("find-tel", cmd_find_tel, "Motif scan + coordinates (no plots)"),
        ("coverage", cmd_coverage, "ONT/HiFi telomere bin depth"),
        ("alignments", cmd_alignments, "Spanning reads from BAM"),
        ("plot", cmd_plot, "All PDFs in one pass"),
        ("all", cmd_all, "Compute then plot"),
    )
    for name, func, help_txt in steps:
        sp = sub.add_parser(name, help=help_txt)
        _add_shared_args(sp)
        sp.set_defaults(func=func)
    return p


def main() -> int:
    args = build_main_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
