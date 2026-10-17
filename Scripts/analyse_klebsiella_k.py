#!/usr/bin/env python3
"""
Summarise Kaptive 3.x Klebsiella K-locus results.
Match confidence is binary in Kaptive 3.x: 'Typeable' or 'Untypeable'.
"""

import csv
from collections import Counter, defaultdict

FILE = "Kaptive_K_results_with_ST.tsv"

# Column indices (0-based)
C_ASSEMBLY   = 3
C_ST         = 4
C_LOCUS      = 5
C_TYPE       = 6    # "Best match type" (e.g. "Capsule null")
C_CONFIDENCE = 7    # "Typeable" / "Untypeable"
C_IDENTITY   = 9
C_COVERAGE   = 10


def load(path):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for r in reader:
            if len(r) < 11:
                continue
            rows.append({
                "assembly":   r[C_ASSEMBLY].strip(),
                "st":         r[C_ST].strip(),
                "locus":      r[C_LOCUS].strip(),
                "match_type": r[C_TYPE].strip(),
                "confidence": r[C_CONFIDENCE].strip(),
                "identity":   r[C_IDENTITY].strip(),
                "coverage":   r[C_COVERAGE].strip(),
            })
    return rows


def main():
    rows = load(FILE)
    print(f"File: {FILE}")
    print(f"Total rows: {len(rows)}")
    print()

    # ---- confidence breakdown ----
    print("Match confidence values present in the file:")
    for k, v in Counter(r["confidence"] for r in rows).most_common():
        print(f"  {k:<15} {v:>5}")
    print()

    # ---- typeability per ST ----
    sts = sorted({r["st"] for r in rows})
    print("Typeability per ST:")
    print(f"  {'ST':<8}{'total':>7}{'typeable':>10}{'%typeable':>12}")
    for st in sts:
        sub = [r for r in rows if r["st"] == st]
        t = sum(1 for r in sub if r["confidence"] == "Typeable")
        pct = 100.0 * t / len(sub) if sub else 0
        print(f"  {st:<8}{len(sub):>7}{t:>10}{pct:>11.1f}%")
    print()

    # ---- locus counts per ST, only Typeable rows ----
    by_st = defaultdict(Counter)
    for r in rows:
        if r["confidence"] == "Typeable":
            by_st[r["st"]][r["locus"]] += 1

    print("Locus counts per ST (Typeable rows only):")
    for st in sts:
        c = by_st.get(st)
        if not c:
            print(f"\n  {st}: no Typeable loci")
            continue
        total = sum(c.values())
        print(f"\n  {st}  ({total} typeable genomes)")
        for locus, k in c.most_common():
            pct = 100.0 * k / total
            print(f"    {locus:<12} {k:>4}  ({pct:>5.1f}%)")

    # ---- global top loci ----
    print()
    print("Top loci across all STs (Typeable only):")
    total = Counter()
    for st in by_st:
        total.update(by_st[st])
    for locus, k in total.most_common(30):
        print(f"  {locus:<12} {k:>4}")


if __name__ == "__main__":
    main()