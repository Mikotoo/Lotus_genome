#!/usr/bin/env python3
"""
Parse barrnap rRNA GFF files to count 18S/5.8S/28S/5S features per rDNA array,
and compute the number of complete 45S repeat units.

Output: rDNA array composition table for Gifu and MG20.
"""

import sys
from collections import defaultdict
from pathlib import Path

# ── rDNA array definitions (1-based coordinates, matching GFF) ──────────────
# Format: label -> (chrom, start, end, expected_type)
# expected_type: "45S" or "5S"

GIFU_ARRAYS = [
    ("Array-1  (Chr2 45S)",   "Chr2", 537436,   949316,   "45S"),
    ("Array-2  (Chr2 5S)",    "Chr2", 24476203,  25027540, "5S"),
    ("Array-3  (Chr5 45S)",   "Chr5", 14943605,  15344676, "45S"),
    ("Array-4  (Chr6 45S)",   "Chr6", 8844951,   8943185,  "45S"),
    ("Array-5  (Chr6 45S)",   "Chr6", 26125599,  27223832, "45S"),
]

MG20_ARRAYS = [
    ("Array-1  (Chr2 45S)",   "Chr2", 1528120,   2857257,  "45S"),
    ("Array-2  (Chr2 5S)",    "Chr2", 26673364,  27897859, "5S"),
    ("Array-3  (Chr5 45S)",   "Chr5", 15013665,  15659889, "45S"),
    ("Array-4  (Chr6 45S)",   "Chr6", 15057397,  15451961, "45S"),
    ("Array-5  (Chr6 45S)",   "Chr6", 26255620,  26622326, "45S"),
]

# ── GFF paths ───────────────────────────────────────────────────────────────
GIFU_GFF = Path(
    "{PROJ_04_LOTUS_GENOME}/"
    "03_annotation/Gifu/rDNA/Gifu_rRNA.gff"
)
MG20_GFF = Path(
    "{PROJ_04_LOTUS_GENOME}/"
    "03_annotation/MG20/rDNA/MG20_rRNA.gff"
)


def parse_rrna_type(name_attr: str) -> str:
    """Extract rRNA subunit type from barrnap GFF Name attribute."""
    if "5S" in name_attr and "5_8S" not in name_attr and "5.8S" not in name_attr:
        return "5S"
    if "5_8S" in name_attr or "5.8S" in name_attr:
        return "5.8S"
    if "18S" in name_attr:
        return "18S"
    if "28S" in name_attr:
        return "28S"
    return "other"


def load_gff(gff_path: Path):
    """Load barrnap GFF, return list of (chrom, start, end, rrna_type, partial)."""
    features = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, _, ftype, start, end, _, _, _, attr = parts
            if ftype != "rRNA":
                continue
            rrna_type = parse_rrna_type(attr)
            is_partial = "partial" in attr
            features.append((chrom, int(start), int(end), rrna_type, is_partial))
    return features


def count_in_region(features, chrom, region_start, region_end):
    """Count 18S/5.8S/28S/5S features with any overlap to the region."""
    counts = {"18S": 0, "5.8S": 0, "28S": 0, "5S": 0}
    partials = {"18S": 0, "5.8S": 0, "28S": 0, "5S": 0}
    for chrom_f, start, end, rrna_type, is_partial in features:
        if chrom_f != chrom:
            continue
        if start > region_end or end < region_start:
            continue
        if rrna_type in counts:
            counts[rrna_type] += 1
            if is_partial:
                partials[rrna_type] += 1
    return counts, partials


def fmt_complete(total, partial):
    """Format as 'total (N partial)' if partials exist, else just total."""
    if partial > 0:
        return f"{total} ({partial} partial)"
    return str(total)


def process_species(species_name, arrays, gff_path):
    """Print a table for one species."""
    features = load_gff(gff_path)
    print(f"## {species_name}")
    print()
    header = (
        f"{'Array':<28s} {'Chrom':>5s} {'Coord (bp)':>24s} "
        f"{'Length (bp)':>12s} {'18S':>8s} {'5.8S':>8s} "
        f"{'28S':>8s} {'5S':>8s} {'45S units':>10s}"
    )
    print(header)
    print("-" * len(header))

    for label, chrom, rs, re, expected in arrays:
        length = re - rs
        coord_str = f"{rs:,}-{re:,}"
        counts, partials = count_in_region(features, chrom, rs, re)

        if expected == "45S":
            # Complete 45S repeat units = min(18S, 5.8S, 28S)
            units_45s = min(counts["18S"], counts["5.8S"], counts["28S"])
            units_str = str(units_45s)
        else:
            units_str = "—"

        row = (
            f"{label:<28s} {chrom:>5s} {coord_str:>24s} "
            f"{length:>12,d} "
            f"{fmt_complete(counts['18S'], partials['18S']):>8s} "
            f"{fmt_complete(counts['5.8S'], partials['5.8S']):>8s} "
            f"{fmt_complete(counts['28S'], partials['28S']):>8s} "
            f"{fmt_complete(counts['5S'], partials['5S']):>8s} "
            f"{units_str:>10s}"
        )
        print(row)
    print()


def main():
    print("# rDNA array composition table for Lotus japonicus")
    print()
    process_species("Gifu", GIFU_ARRAYS, GIFU_GFF)
    process_species("MG20", MG20_ARRAYS, MG20_GFF)

    # ── TSV output for spreadsheet ──
    tsv_path = Path(
        "{PROJ_04_LOTUS_GENOME}/"
        "03_annotation/rDNA_array_summary.tsv"
    )
    with open(tsv_path, "w") as f:
        hdr = [
            "Species", "Array", "Chromosome", "Start", "End",
            "Length_bp", "18S_count", "18S_partial",
            "5.8S_count", "5.8S_partial",
            "28S_count", "28S_partial",
            "5S_count", "5S_partial",
            "Complete_45S_units",
        ]
        f.write("\t".join(hdr) + "\n")

        for species, arrays, gff in [
            ("Gifu", GIFU_ARRAYS, GIFU_GFF),
            ("MG20", MG20_ARRAYS, MG20_GFF),
        ]:
            feats = load_gff(gff)
            for label, chrom, rs, re, expected in arrays:
                counts, partials = count_in_region(feats, chrom, rs, re)
                if expected == "45S":
                    units = min(counts["18S"], counts["5.8S"], counts["28S"])
                else:
                    units = 0
                f.write("\t".join(map(str, [
                    species,
                    label.strip().split("(")[0].strip(),
                    chrom, rs, re,
                    re - rs,
                    counts["18S"], partials["18S"],
                    counts["5.8S"], partials["5.8S"],
                    counts["28S"], partials["28S"],
                    counts["5S"], partials["5S"],
                    units,
                ])) + "\n")

    print(f"TSV saved to: {tsv_path}")


if __name__ == "__main__":
    main()
