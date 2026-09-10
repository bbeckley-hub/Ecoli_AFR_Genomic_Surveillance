import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

mut = pd.read_csv('mutations.csv')
samples = pd.read_csv('sample_overview.csv')

clinically_relevant = ['gyrA_S83L', 'gyrA_D87N', 'gyrA_S83A', 'gyrA_D87G',
                       'parC_S80I', 'parC_E84G', 'parC_E84V', 'parC_E84K',
                       'parE_I529L', 'parE_S458A',
                       'pmrB_D283G', 'pmrB_E123D', 'pmrB_H2R', 'pmrB_V351I',
                       'uhpT_E350Q']

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

presence['ST'] = presence.index.map(samples.set_index('Sample')['ST'].to_dict())
presence = presence.dropna(subset=['ST'])
presence['ST'] = presence['ST'].astype(str)

st_list = ['ST10', 'ST131', 'ST167']
results = []
for key in gene_genomes.keys():
    counts = {}
    for st in st_list:
        sub = presence[presence['ST'] == st]
        present = sub[key].sum()
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
            'Mutation': key,
            'ST': st,
            'Present_ST': table[i][0],
            'Total_ST': table[i][0] + table[i][1],
            'OddsRatio': or_val,
            'p_fisher': p_fish
        })

df_res = pd.DataFrame(results)
df_res['p_adj'] = multipletests(df_res['p_fisher'], method='fdr_bh')[1]
df_res_sig = df_res[df_res['p_adj'] < 0.05].sort_values('p_adj')

print("Significant clinical mutation-ST associations (adj p < 0.05):")
print(df_res_sig[['Mutation', 'ST', 'Present_ST', 'Total_ST', 'OddsRatio', 'p_adj']].to_string(index=False))
