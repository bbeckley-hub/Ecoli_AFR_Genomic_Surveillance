import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests


# ============================================================
# 1. LOAD DATA
# ============================================================

INPUT = "ecoli_comprehensive.csv"

df = pd.read_csv(INPUT, sep="\t")

print("=" * 70)
print("DATASET")
print("=" * 70)
print(f"Number of genomes: {len(df)}")


# ============================================================
# 2. HARMONIZE COUNTRY NAMES
# ============================================================

df["Country"] = (
    df["Country"]
    .astype(str)
    .str.strip()
    .replace({
        "Ethopia": "Ethiopia",
        "kenya": "Kenya"
    })
)

print("\nCountry labels after harmonization:")
print(df["Country"].value_counts().sort_index())


# ============================================================
# 3. DESCRIPTIVE STATISTICS
# ============================================================

categorical_vars = [
    "Country",
    "MLST",
    "Phylogroup",
    "Serotype",
    "CH_Type"
]

for var in categorical_vars:
    print("\n" + "=" * 70)
    print(var)
    print("=" * 70)
    print(df[var].value_counts(dropna=False))


# ============================================================
# 4. FUNCTION FOR CHI-SQUARE + CRAMER'S V
# ============================================================

def chi_square_test(df, var1, var2):

    table = pd.crosstab(df[var1], df[var2])

    chi2, p, dof, expected = chi2_contingency(table)

    n = table.values.sum()
    r, k = table.shape

    # Cramer's V
    phi2 = chi2 / n

    # Bias-corrected Cramer's V
    phi2corr = max(
        0,
        phi2 - ((k - 1) * (r - 1)) / (n - 1)
    )

    rcorr = r - ((r - 1) ** 2) / (n - 1)
    kcorr = k - ((k - 1) ** 2) / (n - 1)

    cramers_v = np.sqrt(
        phi2corr / min((kcorr - 1), (rcorr - 1))
    )

    print("\n" + "=" * 70)
    print(f"{var1} × {var2}")
    print("=" * 70)

    print("\nContingency table:")
    print(table)

    print("\nChi-square:")
    print(f"χ² = {chi2:.4f}")

    print(f"df = {dof}")

    print(f"P = {p:.6g}")

    print(f"Cramer's V = {cramers_v:.4f}")

    print("\nExpected counts:")
    print(
        pd.DataFrame(
            expected,
            index=table.index,
            columns=table.columns
        )
    )

    return table, chi2, dof, p, cramers_v


# ============================================================
# 5. ST × COUNTRY
# ============================================================

chi_square_test(
    df,
    "MLST",
    "Country"
)


# ============================================================
# 6. ST × PHYLOGROUP
# ============================================================

chi_square_test(
    df,
    "MLST",
    "Phylogroup"
)


# ============================================================
# 7. ST × SEROTYPE
# ============================================================

chi_square_test(
    df,
    "MLST",
    "Serotype"
)


# ============================================================
# 8. ST × CH TYPE
# ============================================================

chi_square_test(
    df,
    "MLST",
    "CH_Type"
)


# ============================================================
# 9. STANDARDIZED RESIDUALS
# ============================================================
#
# These tell us WHICH cells drive a significant association.
#
# Values approximately:
# > +2   over-represented
# < -2   under-represented
#
# ============================================================

def standardized_residuals(df, var1, var2):

    table = pd.crosstab(df[var1], df[var2])

    chi2, p, dof, expected = chi2_contingency(table)

    residuals = (
        table.values - expected
    ) / np.sqrt(expected)

    residual_df = pd.DataFrame(
        residuals,
        index=table.index,
        columns=table.columns
    )

    print("\n" + "=" * 70)
    print(f"STANDARDIZED RESIDUALS: {var1} × {var2}")
    print("=" * 70)

    print(residual_df.round(2))

    return residual_df


standardized_residuals(df, "MLST", "Country")
standardized_residuals(df, "MLST", "Phylogroup")
standardized_residuals(df, "MLST", "Serotype")
standardized_residuals(df, "MLST", "CH_Type")


# ============================================================
# 10. ASSEMBLY METRICS
# ============================================================

print("\n" + "=" * 70)
print("ASSEMBLY METRICS")
print("=" * 70)

metrics = [
    "N50",
    "GC_Content",
    "Total_Length"
]

print(df[metrics].describe())


print("\nAssembly metrics by ST:")
print(
    df.groupby("MLST")[metrics]
    .agg(["count", "median", "mean", "std", "min", "max"])
)


# ============================================================
# 11. SAVE CLEANED METADATA
# ============================================================

df.to_csv(
    "ecoli_comprehensive_harmonized.tsv",
    sep="\t",
    index=False
)

print("\n")
print("=" * 70)
print("Analysis complete.")
print("Harmonized dataset saved as:")
print("ecoli_comprehensive_harmonized.tsv")
print("=" * 70)
