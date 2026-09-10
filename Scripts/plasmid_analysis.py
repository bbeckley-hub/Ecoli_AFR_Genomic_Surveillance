import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

plas = pd.read_csv('plasmid_markers.csv')
samples = pd.read_csv('sample_overview.csv')

plas = plas[plas['Database'] == 'PLASMIDFINDER']

mapping = {
    'Col156_1': 'Col156',
    'Col(pHAD28)_1': 'ColpHAD28',
    'Col(BS512)_1': 'ColBS512',
    'Col(MG828)_1': 'ColMG828',
    'ColRNAI_1': 'ColRNAI',
    'ColpVC_1': 'ColpVC',
    'IncFIA_1': 'IncF',
    'IncFIB(AP001918)_1': 'IncF',
    'IncFIC(FII)_1': 'IncF',
    'IncFII_1': 'IncF',
    'IncFII(pRSB107)_1_pRSB107': 'IncF',
    'IncI1_1_Alpha': 'IncI1',
    'IncX3_1': 'IncX3',
    'IncY_1': 'IncY',
    'IncN4_1': 'IncN4',
    'IncFII(29)_1_pUTI89': 'IncF',
    'IncFIB(H89-PhagePlasmid)_1': 'IncF',
    'IncFIB_1': 'IncF',
    'IncFIC_1': 'IncF',
}

for raw in plas['Marker'].unique():
    if raw not in mapping:
        if raw.startswith('Inc') or raw.startswith('Col'):
            mapping[raw] = raw
        else:
            mapping[raw] = raw

plas['Canonical'] = plas['Marker'].map(mapping)

gene_genomes = {}
for _, row in plas.iterrows():
    canon = row['Canonical']
    genomes_str = row['Genomes']
    if pd.isna(genomes_str):
        continue
    genomes = set(genomes_str.split(';'))
    if canon not in gene_genomes:
        gene_genomes[canon] = set()
    gene_genomes[canon].update(genomes)

canonical_replicons = list(gene_genomes.keys())
print(f"Canonical replicon types: {len(canonical_replicons)}")
print(sorted(canonical_replicons))

all_genomes = samples['Sample'].tolist()
presence = pd.DataFrame(0, index=all_genomes, columns=canonical_replicons)

for canon, genomes_set in gene_genomes.items():
    for g in genomes_set:
        if g in presence.index:
            presence.loc[g, canon] = 1

presence['ST'] = presence.index.map(samples.set_index('Sample')['ST'].to_dict())
presence = presence.dropna(subset=['ST'])
presence['ST'] = presence['ST'].astype(str)

st_list = ['ST10', 'ST131', 'ST167']

print("\nPlasmid replicon prevalence by ST:")
for canon in canonical_replicons:
    counts = []
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[canon].sum()
        total = len(sub)
        counts.append(f"{present}/{total} ({100*present/total:.1f}%)")
    print(f"{canon}: ST10 {counts[0]}, ST131 {counts[1]}, ST167 {counts[2]}")

plas_burden = presence[canonical_replicons].sum(axis=1)
samples['Plas_Count'] = samples['Sample'].map(plas_burden.to_dict())

def q1(x): return x.quantile(0.25)
def q3(x): return x.quantile(0.75)

st_summary = samples.groupby('ST')['Plas_Count'].agg(
    median='median', Q1=q1, Q3=q3, N='count'
)
print("\nPlasmid burden per ST:")
print(st_summary)

groups = [samples[samples['ST'] == st]['Plas_Count'].values for st in st_list]
h_stat, p_kw = kruskal(*groups)
print(f"\nKruskal-Wallis: H = {h_stat:.2f}, p = {p_kw:.4e}")

pairwise = []
for i in range(3):
    for j in range(i+1, 3):
        u_stat, p_mw = mannwhitneyu(groups[i], groups[j], alternative='two-sided')
        pairwise.append((st_list[i], st_list[j], p_mw))

_, p_adj, _, _ = multipletests([p for _, _, p in pairwise], method='bonferroni')
print("\nPairwise comparisons (Bonferroni):")
for (st1, st2, _), p_adj_val in zip(pairwise, p_adj):
    print(f"{st1} vs {st2}: adjusted p = {p_adj_val:.4e}")

results = []
for canon in canonical_replicons:
    counts = {}
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[canon].sum()
        absent = len(sub) - present
        counts[st] = (present, absent)
    table = np.array([[counts[st][0], counts[st][1]] for st in st_list])
    for i, st in enumerate(st_list):
        other_present = sum(table[j][0] for j in range(3) if j != i)
        other_absent = sum(table[j][1] for j in range(3) if j != i)
        tbl2 = np.array([[table[i][0], table[i][1]],
                         [other_present, other_absent]])
        or_val, p_fish = fisher_exact(tbl2, alternative='two-sided')
        results.append({
            'Replicon': canon,
            'ST': st,
            'Present_ST': table[i][0],
            'Total_ST': table[i][0] + table[i][1],
            'Present_others': other_present,
            'Total_others': other_present + other_absent,
            'OddsRatio': or_val,
            'p_fisher': p_fish
        })

results_df = pd.DataFrame(results)
results_df['p_fisher_adj'] = multipletests(results_df['p_fisher'], method='fdr_bh')[1]
results_df_sorted = results_df.sort_values('p_fisher_adj')
results_df_sorted.to_csv('plasmid_associations.csv', index=False)

print("\nSignificant plasmid associations (adj p < 0.05):")
sig = results_df_sorted[results_df_sorted['p_fisher_adj'] < 0.05]
print(sig[['Replicon', 'ST', 'Present_ST', 'Total_ST', 'OddsRatio', 'p_fisher_adj']].to_string(index=False))
