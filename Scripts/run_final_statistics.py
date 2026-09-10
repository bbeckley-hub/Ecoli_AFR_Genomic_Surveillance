#!/usr/bin/env python3
"""
Complete statistical analysis for E. coli genome paper.
Handles sparse tables with Monte Carlo simulation.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 1. LOAD DATA
# ============================================================================

# Load comprehensive data
df = pd.read_csv('ecoli_comprehensive.csv', sep='\t')

# Clean sample names (remove .fna if present)
df['Assembly'] = df['Assembly'].str.replace('.fna', '')

print(f"Loaded {len(df)} genomes")

# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

def monte_carlo_chi2(observed, n_simulations=10000, seed=42):
    """
    Monte Carlo simulation for chi-square test on sparse tables.
    Returns: chi2_stat, p_value_monte_carlo, p_value_asymptotic
    """
    np.random.seed(seed)
    
    # Expected counts under independence
    row_sums = observed.sum(axis=1)
    col_sums = observed.sum(axis=0)
    total = observed.sum()
    expected = np.outer(row_sums, col_sums) / total
    
    # Asymptotic chi-square
    chi2_asymp, p_asymp, _, _ = chi2_contingency(observed, correction=False)
    
    # Monte Carlo simulation
    chi2_sim = []
    n_rows, n_cols = observed.shape
    
    for _ in range(n_simulations):
        # Generate random table with same row/col sums
        sim_table = np.zeros_like(observed, dtype=float)
        for i in range(n_rows):
            # Multinomial draw for each row
            probs = col_sums / total
            draw = np.random.multinomial(row_sums[i], probs)
            sim_table[i, :] = draw
        
        # Compute chi-square for this simulation
        sim_chi2 = ((sim_table - expected)**2 / expected).sum()
        chi2_sim.append(sim_chi2)
    
    # Monte Carlo p-value
    p_mc = (np.array(chi2_sim) >= chi2_asymp).mean()
    
    return chi2_asymp, p_asymp, p_mc, expected

def table_assumptions(observed):
    """Check contingency table assumptions."""
    row_sums = observed.sum(axis=1)
    col_sums = observed.sum(axis=0)
    total = observed.sum()
    expected = np.outer(row_sums, col_sums) / total
    
    n_cells = observed.size
    n_expected_lt5 = (expected < 5).sum()
    pct_lt5 = n_expected_lt5 / n_cells * 100
    min_expected = expected.min()
    
    return {
        'n_cells': n_cells,
        'n_expected_lt5': n_expected_lt5,
        'pct_expected_lt5': pct_lt5,
        'min_expected': min_expected,
        'expected': expected
    }

def cramers_v(observed):
    """Calculate Cramér's V effect size."""
    chi2, _, _, _ = chi2_contingency(observed, correction=False)
    n = observed.sum()
    min_dim = min(observed.shape) - 1
    if min_dim == 0:
        return 0
    return np.sqrt(chi2 / (n * min_dim))

def standardize_residuals(observed, expected):
    """Calculate standardized residuals."""
    return (observed - expected) / np.sqrt(expected * (1 - observed.sum(axis=1).reshape(-1,1)/observed.sum()) * (1 - observed.sum(axis=0).reshape(1,-1)/observed.sum()))

def analyze_contingency(observed, row_labels, col_labels, table_name, use_mc=True):
    """Complete analysis for a contingency table."""
    # Assumptions
    assump = table_assumptions(observed)
    
    # Chi-square with Monte Carlo if sparse
    if use_mc and assump['pct_expected_lt5'] > 20:
        chi2, p_asymp, p_mc, expected = monte_carlo_chi2(observed)
        p_final = p_mc
        method = "Monte Carlo (10,000 simulations)"
    else:
        chi2, p_asymp, p_mc, expected = monte_carlo_chi2(observed, n_simulations=100)
        p_final = p_asymp
        method = "Asymptotic chi-square"
    
    # Effect size
    v = cramers_v(observed)
    
    # Residuals
    residuals = standardize_residuals(observed, expected)
    
    # Build results
    results = {
        'table_name': table_name,
        'chi2': chi2,
        'df': (observed.shape[0]-1) * (observed.shape[1]-1),
        'p_value': p_final,
        'p_asymptotic': p_asymp,
        'p_mc': p_mc,
        'method': method,
        'cramers_v': v,
        'n_cells': assump['n_cells'],
        'n_expected_lt5': assump['n_expected_lt5'],
        'pct_expected_lt5': assump['pct_expected_lt5'],
        'min_expected': assump['min_expected'],
        'observed': observed,
        'expected': expected,
        'residuals': residuals,
        'row_labels': row_labels,
        'col_labels': col_labels
    }
    
    return results

# ============================================================================
# 3. RUN ANALYSES
# ============================================================================

results = {}

# 3.1 ST × Country
ct = pd.crosstab(df['MLST'], df['Country'])
row_labels = ct.index.tolist()
col_labels = ct.columns.tolist()
observed = ct.values
results['st_country'] = analyze_contingency(observed, row_labels, col_labels, 'ST × Country')

# 3.2 ST × Phylogroup
ct = pd.crosstab(df['MLST'], df['Phylogroup'])
row_labels = ct.index.tolist()
col_labels = ct.columns.tolist()
observed = ct.values
results['st_phylogroup'] = analyze_contingency(observed, row_labels, col_labels, 'ST × Phylogroup')

# 3.3 ST × Serotype (only include serotypes with >5 occurrences to reduce sparsity)
sero_counts = df['Serotype'].value_counts()
sero_keep = sero_counts[sero_counts >= 5].index.tolist()
df_sero = df[df['Serotype'].isin(sero_keep)]
ct = pd.crosstab(df_sero['MLST'], df_sero['Serotype'])
row_labels = ct.index.tolist()
col_labels = ct.columns.tolist()
observed = ct.values
results['st_serotype'] = analyze_contingency(observed, row_labels, col_labels, 'ST × Serotype (≥5 occurrences)')

# 3.4 ST × CH Type (only include CH types with >5 occurrences)
ch_counts = df['CH_Type'].value_counts()
ch_keep = ch_counts[ch_counts >= 5].index.tolist()
df_ch = df[df['CH_Type'].isin(ch_keep)]
ct = pd.crosstab(df_ch['MLST'], df_ch['CH_Type'])
row_labels = ct.index.tolist()
col_labels = ct.columns.tolist()
observed = ct.values
results['st_chtype'] = analyze_contingency(observed, row_labels, col_labels, 'ST × CH Type (≥5 occurrences)')

# 3.5 ST × Pathotype
#ct = pd.crosstab(df['MLST'], df['Pathotype'])
#row_labels = ct.index.tolist()
#col_labels = ct.columns.tolist()
#observed = ct.values
#results['st_pathotype'] = analyze_contingency(observed, row_labels, col_labels, 'ST × Pathotype')

# 3.6 ST × O-Type (only O-types with >5 occurrences)
o_counts = df['O_Type'].value_counts()
o_keep = o_counts[o_counts >= 5].index.tolist()
df_o = df[df['O_Type'].isin(o_keep)]
ct = pd.crosstab(df_o['MLST'], df_o['O_Type'])
row_labels = ct.index.tolist()
col_labels = ct.columns.tolist()
observed = ct.values
results['st_otype'] = analyze_contingency(observed, row_labels, col_labels, 'ST × O-Type (≥5 occurrences)')

# ============================================================================
# 4. OUTPUT RESULTS
# ============================================================================

print("\n" + "="*80)
print("STATISTICAL ANALYSIS RESULTS")
print("="*80)

summary_data = []

for key, res in results.items():
    print(f"\n{res['table_name']}:")
    print(f"  χ² = {res['chi2']:.4f}, df = {res['df']}, p = {res['p_value']:.2e}")
    print(f"  Method: {res['method']}")
    print(f"  Cramér's V = {res['cramers_v']:.4f}")
    print(f"  Cells with expected <5: {res['n_expected_lt5']}/{res['n_cells']} ({res['pct_expected_lt5']:.1f}%)")
    print(f"  Min expected count: {res['min_expected']:.3f}")
    
    summary_data.append({
        'Comparison': res['table_name'],
        'χ²': f"{res['chi2']:.2f}",
        'df': res['df'],
        'p-value': f"{res['p_value']:.2e}",
        'Method': res['method'],
        "Cramér's V": f"{res['cramers_v']:.4f}",
        '% Expected <5': f"{res['pct_expected_lt5']:.1f}%"
    })

# ============================================================================
# 5. SAVE OUTPUTS
# ============================================================================

# Summary table
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('statistical_summary.tsv', sep='\t', index=False)
print("\nSummary saved to: statistical_summary.tsv")

# Detailed results with residuals for each comparison
for key, res in results.items():
    # Residuals table
    res_df = pd.DataFrame(res['residuals'], 
                          index=res['row_labels'], 
                          columns=res['col_labels'])
    res_df.to_csv(f'{key}_residuals.tsv', sep='\t')
    print(f"Residuals saved to: {key}_residuals.tsv")
    
    # Observed table
    obs_df = pd.DataFrame(res['observed'], 
                          index=res['row_labels'], 
                          columns=res['col_labels'])
    obs_df.to_csv(f'{key}_observed.tsv', sep='\t')
    
    # Expected table
    exp_df = pd.DataFrame(res['expected'], 
                          index=res['row_labels'], 
                          columns=res['col_labels'])
    exp_df.to_csv(f'{key}_expected.tsv', sep='\t')

# ============================================================================
# 6. IDENTIFY IMPORTANT RESIDUALS
# ============================================================================

print("\n" + "="*80)
print("SIGNIFICANT RESIDUALS (>|2|)")
print("="*80)

for key, res in results.items():
    if key != 'st_country':
        continue
    
    print(f"\n{res['table_name']}:")
    resid = res['residuals']
    for i, row in enumerate(res['row_labels']):
        for j, col in enumerate(res['col_labels']):
            if abs(resid[i,j]) > 2:
                direction = "↑" if resid[i,j] > 0 else "↓"
                print(f"  {row} × {col}: {direction} {resid[i,j]:.2f}")

print("\n" + "="*80)
print("DONE!")
print("="*80)
