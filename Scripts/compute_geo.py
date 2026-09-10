
import pandas as pd

f = "ecoli_comprehensive.csv"
df = pd.read_csv(f, sep="\t")

print("\nSHAPE")
print(df.shape)

for col in ["Country", "MLST", "Phylogroup", "Serotype", "CH_Type"]:
    print(f"\n===== {col} =====")
    print(df[col].value_counts(dropna=False).to_string())

print("\n===== ST × COUNTRY =====")
print(pd.crosstab(df["MLST"], df["Country"]).to_string())

print("\n===== ST × PHYLOGROUP =====")
print(pd.crosstab(df["MLST"], df["Phylogroup"]).to_string())

print("\n===== ST × SEROTYPE =====")
print(pd.crosstab(df["MLST"], df["Serotype"]).to_string())

print("\n===== ST × CH TYPE =====")
print(pd.crosstab(df["MLST"], df["CH_Type"]).to_string())

print("\n===== ASSEMBLY METRICS =====")
print(df[["N50", "GC_Content", "Total_Length"]].describe().to_string())

print("\n===== ASSEMBLY METRICS BY ST =====")
print(df.groupby("MLST")[["N50", "GC_Content", "Total_Length"]]
        .describe().to_string())

