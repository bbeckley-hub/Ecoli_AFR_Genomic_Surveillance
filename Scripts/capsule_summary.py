import pandas as pd
import numpy as np

caps = pd.read_csv('ecoli_results_with_ST.tsv', sep='\t', header=0)

print("Columns:", caps.columns.tolist())

# Clean up column names 
caps.columns = caps.columns.str.strip()

# Check if 'ST' is present 
print("Available columns:", caps.columns.tolist())

# Extract relevant columns
keep_cols = ['Assembly', 'Best match locus', 'Best match type', 'Match confidence', 'Identity', 'Coverage', 'ST']
caps_clean = caps[keep_cols].copy()

# Remove "%" from Identity and Coverage and convert to float
caps_clean['Identity'] = caps_clean['Identity'].str.replace('%', '').astype(float)
caps_clean['Coverage'] = caps_clean['Coverage'].str.replace('%', '').astype(float)

# Filter only ST10, ST131, ST167 (to be safe)
caps_clean = caps_clean[caps_clean['ST'].isin(['ST10', 'ST131', 'ST167'])]

# For Typeable ones, keep only those with high confidence (>90% identity and coverage)

print("\nSummary by ST:")
summary = caps_clean.groupby('ST').agg(
    Total=('ST', 'count'),
    Typeable=('Match confidence', lambda x: (x == 'Typeable').sum()),
    Untypeable=('Match confidence', lambda x: (x == 'Untypeable').sum())
)
print(summary)

print("\nTypeable loci by ST:")
typeable = caps_clean[caps_clean['Match confidence'] == 'Typeable']
loci_summary = typeable.groupby(['ST', 'Best match locus']).size().reset_index(name='Count')
print(loci_summary)

print("\nTypeable types by ST:")
types_summary = typeable.groupby(['ST', 'Best match type']).size().reset_index(name='Count')
print(types_summary)

# For ST167, check the most common high-confidence loci
st167_typeable = typeable[typeable['ST'] == 'ST167']
print("\nST167 Typeable loci with identity/coverage:")
print(st167_typeable[['Best match locus', 'Best match type', 'Identity', 'Coverage']].drop_duplicates().to_string(index=False))
