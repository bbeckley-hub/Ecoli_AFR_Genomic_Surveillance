import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count
import itertools

# Load data
rtab = pd.read_csv('gene_presence_absence.Rtab', sep='\t')
samples = pd.read_csv('sample_overview.csv')

# Transpose and prepare
rtab = rtab.T
rtab.columns = rtab.iloc[0]
rtab = rtab.drop(rtab.index[0])
st_map = samples.set_index('Sample')['ST'].to_dict()
rtab['ST'] = rtab.index.map(st_map)
rtab = rtab.dropna(subset=['ST'])
rtab['ST'] = rtab['ST'].astype(str)
gene_cols = [c for c in rtab.columns if c != 'ST']
presence = rtab[gene_cols].astype(int)

# Filter to shell genes (15-95% prevalence) - already ~3345
n_genomes = len(rtab)
gene_counts = presence.sum(axis=0)
shell_genes = gene_counts[(gene_counts >= 0.15 * n_genomes) & (gene_counts < 0.99 * n_genomes)].index.tolist()
print(f"Testing {len(shell_genes)} shell genes.")

st_list = ['ST10', 'ST131', 'ST167']

def test_gene(gene):
    results = []
    counts = {}
    for st in st_list:
        sub = rtab[rtab['ST'] == st]
        present = sub[gene].sum()
        absent = len(sub) - present
        counts[st] = (present, absent)
    for i, st in enumerate(st_list):
        other_present = 0
        other_absent = 0
        for j, st2 in enumerate(st_list):
            if j != i:
                other_present += counts[st2][0]
                other_absent += counts[st2][1]
        if other_present + other_absent > 0:
            tbl = np.array([[counts[st][0], counts[st][1]],
                            [other_present, other_absent]])
            or_val, p_fish = fisher_exact(tbl, alternative='two-sided')
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
    return results

if __name__ == '__main__':
    with Pool(cpu_count()) as pool:
        all_results = pool.map(test_gene, shell_genes)
    # Flatten results
    flat_results = [item for sublist in all_results for item in sublist]
    df_res = pd.DataFrame(flat_results)
    df_res['p_adj'] = multipletests(df_res['p_fisher'], method='fdr_bh')[1]
    df_res_sig = df_res[df_res['p_adj'] < 0.05].sort_values('p_adj')
    df_res_sig.to_csv('st_enriched_genes_fast.csv', index=False)
    print(f"Significant genes: {len(df_res_sig)}")