#!/usr/bin/env python3
"""
iTOL annotations for Klebsiella-like K-loci in E. coli ST10/ST131/ST167.

Reads: Kaptive_K_results_with_ST.tsv
Writes:
  itol_Klebsiella_K_colorstrip.txt   (single categorical strip)
  itol_Klebsiella_K_binary.txt       (9-column binary heatmap)
"""

import csv
import os
import sys
from collections import Counter

KAPTIVE_FILE = "Kaptive_K_results_with_ST.tsv"
STRIP_OUT    = "itol_Klebsiella_K_colorstrip.txt"
BINARY_OUT   = "itol_Klebsiella_K_binary.txt"
STRIP_WIDTH  = 25

# Column indices (0-based)
C_ASSEMBLY   = 3
C_ST         = 4
C_LOCUS      = 5
C_CONFIDENCE = 7

# Top 9 loci (from the analysis). Order = column / legend order.
TOP_LOCI = [
    "KL9",      # ST167 dominant
    "KL103",    # ST167 dominant
    "KL115",    # ST167
    "KL114",    # ST10
    "KL124",    # ST167 (minor ST10)
    "KL137",    # ST10
    "KL30",     # ST167
    "KL158",    # ST10
    "KL35",     # ST167
]

LOCUS_COLORS = {
    "KL9":   "#e6194b",   # red
    "KL103": "#3cb44b",   # green
    "KL115": "#4363d8",   # blue
    "KL114": "#f58231",   # orange
    "KL124": "#911eb4",   # purple
    "KL137": "#46f0f0",   # cyan
    "KL30":  "#f032e6",   # magenta
    "KL158": "#bcf60c",   # lime
    "KL35":  "#fabebe",   # pink
}
OTHERS_COLOR     = "#d3d3d3"   # light gray
UNTYPEABLE_COLOR = "#000000"  # white


def load_kaptive(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: file not found: {path}")
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)                          # header
        for r in reader:
            if len(r) < 8:
                continue
            rows.append({
                "assembly":   r[C_ASSEMBLY].strip(),
                "st":         r[C_ST].strip(),
                "locus":      r[C_LOCUS].strip(),
                "confidence": r[C_CONFIDENCE].strip(),
            })
    return rows


def categorise(row):
    """Return one of: a TOP_LOCI name, 'Others', or 'Untypeable'."""
    if row["confidence"] != "Typeable":
        return "Untypeable"
    if row["locus"] in TOP_LOCI:
        return row["locus"]
    return "Others"


def write_colorstrip(path, rows):
    # legend: TOP_LOCI in listed order, then Others, then Untypeable
    legend_labels = TOP_LOCI + ["Others", "Untypeable"]
    legend_colors = [LOCUS_COLORS[k] for k in TOP_LOCI] + [OTHERS_COLOR, UNTYPEABLE_COLOR]

    with open(path, "w", encoding="utf-8") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tKlebsiella-like K-locus\n")
        f.write("COLOR\t#000000\n")
        f.write(f"STRIP_WIDTH\t{STRIP_WIDTH}\n")
        f.write("MARGIN\t0\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("SHOW_INTERNAL\t0\n")
        f.write("LEGEND_TITLE\tKlebsiella-like K-locus\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
        f.write("DATA\n")
        for row in rows:
            cat = categorise(row)
            color = LOCUS_COLORS.get(cat, OTHERS_COLOR if cat == "Others" else UNTYPEABLE_COLOR)
            f.write(f"{row['assembly']}\t{color}\t{cat}\n")


def write_binary(path, rows):
    """One column per top locus, 1 = present, 0 = absent."""
    n = len(TOP_LOCI)
    with open(path, "w", encoding="utf-8") as f:
        f.write("DATASET_BINARY\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tKlebsiella-like K-loci\n")
        f.write("COLOR\t#000000\n")
        f.write("FIELD_LABELS\t" + "\t".join(TOP_LOCI) + "\n")
        f.write("FIELD_SHAPES\t" + "\t".join(["1"] * n) + "\n")
        f.write("FIELD_COLORS\t" + "\t".join(["#000000"] * n) + "\n")
        f.write("LEGEND_TITLE\tKlebsiella-like K-loci\n")
        f.write("LEGEND_SHAPES\t1\n")
        f.write("LEGEND_COLORS\t#000000\n")
        f.write("LEGEND_LABELS\t1\n")
        f.write("DATA\n")
        for row in rows:
            cells = []
            for locus in TOP_LOCI:
                hit = (row["confidence"] == "Typeable" and row["locus"] == locus)
                cells.append("1" if hit else "0")
            f.write(row["assembly"] + "\t" + "\t".join(cells) + "\n")


def main():
    rows = load_kaptive(KAPTIVE_FILE)
    print(f"Read {len(rows)} samples from {KAPTIVE_FILE}")

    write_colorstrip(STRIP_OUT, rows)
    write_binary(BINARY_OUT, rows)

    # ---- summary ----
    cats = Counter(categorise(r) for r in rows)
    print(f"\nWrote: {STRIP_OUT}")
    print(f"Wrote: {BINARY_OUT}")
    print()
    print(f"{'Category':<12}{'Count':>7}   {'ST breakdown'}")
    for cat in TOP_LOCI + ["Others", "Untypeable"]:
        n = cats.get(cat, 0)
        if cat in TOP_LOCI:
            # show ST breakdown
            st_counts = Counter(r["st"] for r in rows if categorise(r) == cat)
            st_str = ", ".join(f"{k}={v}" for k, v in st_counts.most_common())
        else:
            st_str = ""
        print(f"{cat:<12}{n:>7}   {st_str}")


if __name__ == "__main__":
    main()
