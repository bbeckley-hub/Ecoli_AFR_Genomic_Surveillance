import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import re

amr = pd.read_csv('amr_genes.csv')
samples = pd.read_csv('sample_overview.csv')

# Use the exact same mapping 
mapping = {
    # Beta-lactamases
    'blaCTX-M-15': 'blaCTX-M-15',
    'blaCTX-M-15_1': 'blaCTX-M-15',
    'CTX-M-15': 'blaCTX-M-15',
    '(Bla)blaCTX-M-15': 'blaCTX-M-15',

    'blaNDM-5': 'blaNDM-5',
    'blaNDM-5_1': 'blaNDM-5',
    'NDM-5': 'blaNDM-5',
    '(Bla)blaNDM-5': 'blaNDM-5',

    'blaOXA-1': 'blaOXA-1',
    'blaOXA-1_1': 'blaOXA-1',
    'OXA-1': 'blaOXA-1',
    '(Bla)blaOXA-1': 'blaOXA-1',

    'blaTEM-1': 'blaTEM-1',
    'blaTEM-1_1': 'blaTEM-1',
    'TEM-1': 'blaTEM-1',
    '(Bla)blaTEM-105': 'blaTEM-1',
    'blaTEM-1B_1': 'blaTEM-1',

    'blaTEM-35': 'blaTEM-35',
    '(Bla)blaTEM-158': 'blaTEM-35',

    'blaCMY-2': 'blaCMY-2',
    'blaCMY-2_1': 'blaCMY-2',
    'CMY-2': 'blaCMY-2',
    '(Bla)blaCMY-111': 'blaCMY-2',

    'blaOXA-48': 'blaOXA-48',
    'blaOXA-48_1': 'blaOXA-48',
    'OXA-48': 'blaOXA-48',
    '(Bla)blaOXA-48': 'blaOXA-48',

    'blaDHA-1': 'blaDHA-1',
    'blaDHA-1_1': 'blaDHA-1',
    '(Bla)blaDHA-1': 'blaDHA-1',

    'blaCTX-M-27': 'blaCTX-M-27',
    'blaCTX-M-27_1': 'blaCTX-M-27',
    '(Bla)blaCTX-M-27': 'blaCTX-M-27',

    'blaCTX-M-55': 'blaCTX-M-55',
    'blaCTX-M-55_1': 'blaCTX-M-55',
    '(Bla)blaCTX-M-55': 'blaCTX-M-55',

    'blaCTX-M-14': 'blaCTX-M-14',
    '(Bla)blaCTX-M-14': 'blaCTX-M-14',

    # Aminoglycosides
    'aac(3)-IId': 'aac(3)-IId',
    'aac(3)-IId_1': 'aac(3)-IId',
    'AAC(3)-IId': 'aac(3)-IId',
    '(AGly)aac3-IId': 'aac(3)-IId',

    "aac(6')-Ib-cr": "aac(6')-Ib-cr",
    "aac(6')-Ib-cr5": "aac(6')-Ib-cr",
    "aac(6')-Ib-cr_1": "aac(6')-Ib-cr",
    "AAC(6')-Ib-cr": "aac(6')-Ib-cr",

    "aph(3'')-Ib": "aph(3'')-Ib",
    "aph(3'')-Ib_5": "aph(3'')-Ib",
    "APH(3'')-Ib": "aph(3'')-Ib",

    "aph(6)-Id": "aph(6)-Id",
    "aph(6)-Id_1": "aph(6)-Id",
    "APH(6)-Id": "aph(6)-Id",

    'aadA1': 'aadA1',
    'aadA1_1': 'aadA1',
    '(AGly)aadA1': 'aadA1',

    'aadA2': 'aadA2',
    'aadA2_1': 'aadA2',
    '(AGly)aadA2': 'aadA2',
    'aadA2_2': 'aadA2',

    'aadA5': 'aadA5',
    'aadA5_1': 'aadA5',
    '(AGly)aadA5': 'aadA5',

    'aac(3)-IIe': 'aac(3)-IIe',
    'AAC(3)-IIe': 'aac(3)-IIe',
    '(AGly)aac3-IIe': 'aac(3)-IIe',

    # Sulphonamides
    'sul1': 'sul1',
    'sul1_5': 'sul1',
    '(Sul)sul1': 'sul1',

    'sul2': 'sul2',
    'sul2_2': 'sul2',
    '(Sul)sul2': 'sul2',

    'sul3': 'sul3',
    '(Sul)sul3': 'sul3',

    # Trimethoprim
    'dfrA1': 'dfrA1',
    'dfrA1_8': 'dfrA1',
    '(Tmt)dfrA1': 'dfrA1',

    'dfrA17': 'dfrA17',
    'dfrA17_1': 'dfrA17',
    '(Tmt)dfrA17': 'dfrA17',

    'dfrA14': 'dfrA14',
    'dfrA14_5': 'dfrA14',
    '(Tmt)dfrA14': 'dfrA14',

    'dfrA12': 'dfrA12',
    'dfrA12_8': 'dfrA12',
    '(Tmt)dfrA12': 'dfrA12',

    'dfrA8': 'dfrA8',
    '(Tmt)dfrA8': 'dfrA8',

    'dfrA7': 'dfrA7',
    '(Tmt)dfrA7': 'dfrA7',

    'dfrA5': 'dfrA5',
    '(Tmt)dfrA5': 'dfrA5',

    # Tetracyclines
    'tet(A)': 'tet(A)',
    'tet(A)_6': 'tet(A)',
    '(Tet)tetA': 'tet(A)',
    'TETA': 'tet(A)',

    'tet(B)': 'tet(B)',
    'tet(B)_2': 'tet(B)',
    '(Tet)tetB': 'tet(B)',
    'TETB': 'tet(B)',

    'tet(M)': 'tet(M)',
    'tet(M)_10': 'tet(M)',
    '(Tet)tetM': 'tet(M)',
    'TETM': 'tet(M)',

    'tet(C)': 'tet(C)',
    'tet(C)_3': 'tet(C)',
    '(Tet)tetC': 'tet(C)',
    'TETC': 'tet(C)',

    'tet(D)': 'tet(D)',
    'tet(D)_1': 'tet(D)',
    '(Tet)tetD': 'tet(D)',
    'TETD': 'tet(D)',

    # Phenicols
    'catB3': 'catB3',
    'catB3_2': 'catB3',
    '(Phe)catB3': 'catB3',

    'catA1': 'catA1',
    'catA1_1': 'catA1',
    '(Phe)catA1': 'catA1',

    'catA2': 'catA2',
    'catA2_1': 'catA2',
    '(Phe)catA2': 'catA2',

    # Macrolides
    'mph(A)': 'mph(A)',
    '(MLS)mph(A)': 'mph(A)',
    'mph(A)_2': 'mph(A)',
    'mph(A)_1': 'mph(A)',
    'MPHA': 'mph(A)',

    'mph(B)': 'mph(B)',
    'mph(B)_1': 'mph(B)',
    '(MLS)mph(B)': 'mph(B)',

    'erm(B)': 'erm(B)',
    'erm(B)_1': 'erm(B)',
    'erm(B)_18': 'erm(B)',
    '(MLS)erm(B)': 'erm(B)',
    'ERMB': 'erm(B)',
    'ErmB': 'erm(B)',

    'rmtB': 'rmtB',
    'rmtB1': 'rmtB',
    'RMTB': 'rmtB',
    '(AGly)rmtB': 'rmtB',

    # Quinolones
    'qnrS1': 'qnrS1',
    'qnrS1_1': 'qnrS1',
    'QnrS1': 'qnrS1',
    'QNRS': 'qnrS1',
    '(Flq)qnr-S1': 'qnrS1',

    'qepA': 'qepA',
    'qepA1': 'qepA',
    'qepA1_1': 'qepA',
    'QepA2': 'qepA',
    '(Flq)qepA': 'qepA',

    'qepA4': 'qepA4',
    'qepA4_1': 'qepA4',

    'qnrB1': 'qnrB1',
    'qnrB4': 'qnrB4',
    'QnrB4': 'qnrB4',
    '(Flq)qnrB1': 'qnrB1',
    '(Flq)qnrB4': 'qnrB4',

    # Colistin
    'mcr-1': 'mcr-1',

    # Others
    'floR': 'floR',
    'floR_4': 'floR',
    '(Phe)floR': 'floR',
    'FLOR': 'floR',

    'fosA': 'fosA',

    'cmlA': 'cmlA',
    'cmlA1': 'cmlA',
    '(Phe)cmlA1': 'cmlA',
}

# Auto-add trailing digit variants
additional = {}
for raw in list(mapping.keys()):
    if re.search(r'_\d+$', raw):
        base = re.sub(r'_\d+$', '', raw)
        if base in mapping:
            additional[raw] = mapping[base]
mapping.update(additional)

# Build gene -> set of genomes (union across databases)
gene_genomes = {}
for _, row in amr.iterrows():
    raw = row['Gene']
    genomes_str = row['Genomes']
    if pd.isna(genomes_str):
        continue
    genomes = set(genomes_str.split(';'))
    canon = mapping.get(raw)
    if canon is None:
        cleaned = re.sub(r'[_\d]+$', '', raw)
        canon = mapping.get(cleaned)
    if canon is not None:
        if canon not in gene_genomes:
            gene_genomes[canon] = set()
        gene_genomes[canon].update(genomes)

# Count per genome (only those in samples)
all_genomes = samples['Sample'].tolist()
amr_count = {g: 0 for g in all_genomes}
for gene, genomes_set in gene_genomes.items():
    for g in genomes_set:
        if g in amr_count:
            amr_count[g] += 1

# Add to samples
samples['AMR_Count'] = samples['Sample'].map(amr_count)

# Compute summary per ST using quantile for Q1 and Q3
def q1(x): return x.quantile(0.25)
def q3(x): return x.quantile(0.75)

st_counts = samples.groupby('ST')['AMR_Count'].agg(
    median='median',
    Q1=q1,
    Q3=q3,
    N='count'
)
print("\nAMR burden per ST:")
print(st_counts)

# Kruskal-Wallis
groups = [samples[samples['ST'] == st]['AMR_Count'].values for st in ['ST10', 'ST131', 'ST167']]
h_stat, p_kw = kruskal(*groups)
print(f"\nKruskal-Wallis: H = {h_stat:.2f}, p = {p_kw:.4e}")

# Dunn's post-hoc (pairwise Mann-Whitney with Bonferroni correction)
st_list = ['ST10', 'ST131', 'ST167']
pairwise = []
for i in range(3):
    for j in range(i+1, 3):
        u_stat, p_mw = mannwhitneyu(groups[i], groups[j], alternative='two-sided')
        pairwise.append((st_list[i], st_list[j], p_mw))
_, p_adj, _, _ = multipletests([p for _, _, p in pairwise], method='bonferroni')
print("\nPairwise comparisons (Bonferroni-adjusted):")
for (st1, st2, _), p_adj_val in zip(pairwise, p_adj):
    print(f"{st1} vs {st2}: adjusted p = {p_adj_val:.4e}")
