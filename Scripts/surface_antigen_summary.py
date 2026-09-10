import pandas as pd
import numpy as np

# Load files
ecoli_k = pd.read_csv('ecoli_results_with_ST.tsv', sep='\t', header=0)
k_k = pd.read_csv('k_results_with_ST.tsv', sep='\t', header=0)
o_k = pd.read_csv('o_results_with_ST.tsv', sep='\t', header=0)

# Clean column names (remove trailing spaces)
for df in [ecoli_k, k_k, o_k]:
    df.columns = df.columns.str.strip()
    # The last column is 'ST' but might be named 'ST.1' if duplicate; we'll rename
    if 'ST.1' in df.columns:
        df['ST'] = df['ST.1']
        df.drop(columns=['ST.1'], inplace=True)

# Keep only the three STs of interest
st_list = ['ST10', 'ST131', 'ST167']
ecoli_k = ecoli_k[ecoli_k['ST'].isin(st_list)]
k_k = k_k[k_k['ST'].isin(st_list)]
o_k = o_k[o_k['ST'].isin(st_list)]

# Helper to summarise
def summarise(df, name):
    print(f"\n{'='*60}")
    print(f"{name} summary")
    print(f"{'='*60}")
    # Typeability
    summary = df.groupby('ST').agg(
        Total=('ST', 'count'),
        Typeable=('Match confidence', lambda x: (x == 'Typeable').sum()),
        Untypeable=('Match confidence', lambda x: (x == 'Untypeable').sum())
    )
    print(summary)
    # For Typeable, top loci
    typeable = df[df['Match confidence'] == 'Typeable']
    if not typeable.empty:
        loci = typeable.groupby(['ST', 'Best match locus']).size().reset_index(name='Count')
        print(f"\nTypeable loci by ST:")
        print(loci)
        types = typeable.groupby(['ST', 'Best match type']).size().reset_index(name='Count')
        print(f"\nTypeable types by ST:")
        print(types)
    else:
        print("\nNo Typeable matches found.")

# Run for each
summarise(ecoli_k, "E. coli Group 2/3 K-loci")
summarise(k_k, "Klebsiella K-loci")
summarise(o_k, "Klebsiella O-loci")

# Specifically check ST167 in Klebsiella K-loci
st167_k = k_k[k_k['ST'] == 'ST167']
print("\n\nST167 Klebsiella K-loci details:")
if not st167_k.empty:
    print(st167_k[['Assembly', 'Best match locus', 'Best match type', 'Match confidence', 'Identity', 'Coverage']].head(10).to_string(index=False))
else:
    print("No ST167 genomes in Klebsiella K-loci file.")

# Check ST167 in Klebsiella O-loci
st167_o = o_k[o_k['ST'] == 'ST167']
print("\nST167 Klebsiella O-loci details (first 10):")
if not st167_o.empty:
    print(st167_o[['Assembly', 'Best match locus', 'Best match type', 'Match confidence', 'Identity', 'Coverage']].head(10).to_string(index=False))
else:
    print("No ST167 genomes in Klebsiella O-loci file.")
