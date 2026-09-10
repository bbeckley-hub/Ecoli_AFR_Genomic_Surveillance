import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

vir = pd.read_csv('virulence_genes.csv')
samples = pd.read_csv('sample_overview.csv')

pathotype_genes = [
    'eae', 'espA', 'espB', 'espD', 'espF', 'espG', 'espJ', 'espV',
    'stx1A', 'stx1B', 'stx2A', 'stx2B',
    'aafA', 'aafD', 'aatA', 'aap',
    'papA', 'papC', 'papG', 'papE', 'papF', 'papH',
    'cnf1', 'hlyA', 'hlyC', 'hlyD',
    'iha', 'fimH',
    'iucA', 'iucB', 'iucC', 'iucD', 'iutA',
    'iroB', 'iroC', 'iroD', 'iroE', 'iroN',
    'fyuA', 'irp2', 'ireA',
    'chuA', 'chuS', 'chuT', 'chuU', 'chuV', 'chuW', 'chuX',
    'kpsU', 'kpsE', 'kpsD', 'kpsT',
    'iss', 'sat', 'pic', 'ibeA',
    'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG',
    'clbH', 'clbI', 'clbJ', 'clbK', 'clbL', 'clbM', 'clbN',
    'clbO', 'clbP', 'clbQ', 'clbS'
]

gene_genomes = {}
for _, row in vir.iterrows():
    gene = row['Gene']
    genomes_str = row['Genomes']
    if pd.isna(genomes_str):
        continue
    genomes = set(genomes_str.split(';'))
    gene_genomes[gene] = genomes

available_genes = [g for g in pathotype_genes if g in gene_genomes]
print(f"Available genes: {len(available_genes)}")

all_genomes = samples['Sample'].tolist()
presence = pd.DataFrame(0, index=all_genomes, columns=available_genes)

for gene in available_genes:
    for g in gene_genomes[gene]:
        if g in presence.index:
            presence.loc[g, gene] = 1

presence['ST'] = presence.index.map(samples.set_index('Sample')['ST'].to_dict())
presence = presence.dropna(subset=['ST'])
presence['ST'] = presence['ST'].astype(str)

st_list = ['ST10', 'ST131', 'ST167']

print("\nVirulence gene prevalence by ST:")
for gene in available_genes:
    counts = []
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[gene].sum()
        total = len(sub)
        counts.append(f"{present}/{total} ({100*present/total:.1f}%)")
    print(f"{gene}: ST10 {counts[0]}, ST131 {counts[1]}, ST167 {counts[2]}")

vir_burden = presence[available_genes].sum(axis=1)
samples['Vir_Count'] = samples['Sample'].map(vir_burden.to_dict())

def q1(x): return x.quantile(0.25)
def q3(x): return x.quantile(0.75)

st_summary = samples.groupby('ST')['Vir_Count'].agg(
    median='median', Q1=q1, Q3=q3, N='count'
)
print("\nVirulence burden per ST:")
print(st_summary)

groups = [samples[samples['ST'] == st]['Vir_Count'].values for st in st_list]
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
for gene in available_genes:
    counts = {}
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[gene].sum()
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
            'Gene': gene,
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
results_df_sorted.to_csv('virulence_associations.csv', index=False)

print("\nSignificant virulence gene associations (adj p < 0.05):")
sig = results_df_sorted[results_df_sorted['p_fisher_adj'] < 0.05]
print(sig[['Gene', 'ST', 'Present_ST', 'Total_ST', 'OddsRatio', 'p_fisher_adj']].to_string(index=False))