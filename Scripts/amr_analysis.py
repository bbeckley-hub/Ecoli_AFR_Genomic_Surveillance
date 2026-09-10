import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import re

amr = pd.read_csv('amr_genes.csv')
samples = pd.read_csv('sample_overview.csv')

acquired_patterns = [
    r'^bla', r'^aac', r'^aph', r'^aad', r'^sul', r'^dfr',
    r'^tet', r'^qnr', r'^qep', r'^mph', r'^cat', r'^mcr',
    r'^floR', r'^fosA', r'^cmlA', r'^erm', r'^mef', r'^msr'
]

def canonical_name(raw):
    raw_clean = re.sub(r'[_\s]+', '', raw)
    raw_lower = raw_clean.lower()
    for pat in acquired_patterns:
        if re.match(pat, raw_lower, re.IGNORECASE):
            base = re.sub(r'_\d+$', '', raw_clean)
            if base.upper() == 'MPHA':
                return 'mph(A)'
            if base.upper() == 'MPHB':
                return 'mph(B)'
            if base.upper() == 'QEPA4':
                return 'qepA4'
            if base.upper() == 'QEPA':
                return 'qepA'
            return base
    return None

mapping = {}
for _, row in amr.iterrows():
    raw = row['Gene']
    canon = canonical_name(raw)
    if canon:
        mapping[raw] = canon

gene_genomes = {}
for _, row in amr.iterrows():
    raw = row['Gene']
    genomes_str = row['Genomes']
    if pd.isna(genomes_str):
        continue
    genomes = set(genomes_str.split(';'))
    canon = mapping.get(raw)
    if canon:
        if canon not in gene_genomes:
            gene_genomes[canon] = set()
        gene_genomes[canon].update(genomes)

canonical_genes = list(gene_genomes.keys())
all_genomes = samples['Sample'].tolist()
presence = pd.DataFrame(0, index=all_genomes, columns=canonical_genes)

for gene, genomes_set in gene_genomes.items():
    for g in genomes_set:
        if g in presence.index:
            presence.loc[g, gene] = 1

st_map = samples.set_index('Sample')['ST'].to_dict()
presence['ST'] = presence.index.map(st_map)
presence = presence.dropna(subset=['ST'])
presence['ST'] = presence['ST'].astype(str)

results = []
st_list = ['ST10', 'ST131', 'ST167']

for gene in canonical_genes:
    counts = {}
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[gene].sum()
        absent = len(sub) - present
        counts[st] = (present, absent)
    for i, st in enumerate(st_list):
        other_present = sum(counts[st2][0] for st2 in st_list if st2 != st)
        other_absent = sum(counts[st2][1] for st2 in st_list if st2 != st)
        tbl2 = np.array([[counts[st][0], counts[st][1]],
                         [other_present, other_absent]])
        or_val, p_fish = fisher_exact(tbl2, alternative='two-sided')
        results.append({
            'Gene': gene,
            'ST': st,
            'Present_ST': counts[st][0],
            'Total_ST': counts[st][0] + counts[st][1],
            'Present_others': other_present,
            'Total_others': other_present + other_absent,
            'OddsRatio': or_val,
            'p_fisher': p_fish
        })

results_df = pd.DataFrame(results)
results_df['p_fisher_adj'] = multipletests(results_df['p_fisher'], method='fdr_bh')[1]
results_df_sorted = results_df.sort_values('p_fisher_adj')
results_df_sorted.to_csv('amr_stat_results.csv', index=False)

print(results_df_sorted[results_df_sorted['p_fisher_adj'] < 0.05][['Gene', 'ST', 'Present_ST', 'Total_ST', 'OddsRatio', 'p_fisher_adj']].to_string(index=False))