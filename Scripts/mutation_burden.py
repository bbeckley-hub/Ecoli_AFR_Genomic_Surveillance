import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

mut = pd.read_csv('mutations.csv')
samples = pd.read_csv('sample_overview.csv')

clinically_relevant_mutations = [
    ('gyrA', 'S83L'), ('gyrA', 'D87N'), ('gyrA', 'S83A'), ('gyrA', 'D87G'),
    ('parC', 'S80I'), ('parC', 'E84G'), ('parC', 'E84V'), ('parC', 'E84K'),
    ('parE', 'I529L'), ('parE', 'S458A'),
    ('pmrB', 'D283G'), ('pmrB', 'E123D'), ('pmrB', 'H2R'), ('pmrB', 'V351I'),
    ('uhpT', 'E350Q')
]

mut['Genomes'] = mut['Genomes'].fillna('')
mut['Genomes_list'] = mut['Genomes'].apply(lambda x: set(x.split(';')) if x else set())

gene_genomes = {}
for _, row in mut.iterrows():
    key = f"{row['Gene']}_{row['Mutation']}"
    gene_genomes[key] = row['Genomes_list']

all_genomes = samples['Sample'].tolist()
presence = pd.DataFrame(0, index=all_genomes, columns=gene_genomes.keys())
for key, genomes_set in gene_genomes.items():
    for g in genomes_set:
        if g in presence.index:
            presence.loc[g, key] = 1

# Subset to clinically relevant ones
relevant_keys = []
for gene, mut_name in clinically_relevant_mutations:
    # Try exact match with "_nan" suffix 
    key1 = f"{gene}_{mut_name}_nan"
    # Also try without _nan (if any)
    key2 = f"{gene}_{mut_name}"
    if key1 in presence.columns:
        relevant_keys.append(key1)
    elif key2 in presence.columns:
        relevant_keys.append(key2)

print(f"Found {len(relevant_keys)} clinically relevant mutations in data")

# Count per genome
presence['Mut_Count'] = presence[relevant_keys].sum(axis=1)

# Add ST
st_map = samples.set_index('Sample')['ST'].to_dict()
presence['ST'] = presence.index.map(st_map)
presence = presence.dropna(subset=['ST'])
presence['ST'] = presence['ST'].astype(str)

# Summarise
def q1(x): return x.quantile(0.25)
def q3(x): return x.quantile(0.75)

st_summary = presence.groupby('ST')['Mut_Count'].agg(
    median='median', Q1=q1, Q3=q3, N='count'
)
print("\nMutation burden per ST:")
print(st_summary)

groups = [presence[presence['ST'] == st]['Mut_Count'].values for st in ['ST10', 'ST131', 'ST167']]
h_stat, p_kw = kruskal(*groups)
print(f"\nKruskal-Wallis: H = {h_stat:.2f}, p = {p_kw:.4e}")

pairwise = []
st_list = ['ST10', 'ST131', 'ST167']
for i in range(3):
    for j in range(i+1, 3):
        u_stat, p_mw = mannwhitneyu(groups[i], groups[j], alternative='two-sided')
        pairwise.append((st_list[i], st_list[j], p_mw))
_, p_adj, _, _ = multipletests([p for _, _, p in pairwise], method='bonferroni')
print("\nPairwise comparisons (Bonferroni):")
for (st1, st2, _), p_adj_val in zip(pairwise, p_adj):
    print(f"{st1} vs {st2}: adjusted p = {p_adj_val:.4e}")
