```python
#!/usr/bin/env python3
"""
Microbiome Data Analysis Pipeline: Taxonomy, Function, and Pathogen Load Interactions
Author: Sergio Castañeda
Description: Comprehensive pipeline for processing harmonized MetaPhlAn (taxonomic) 
             and HUMAnN (functional) profiles, including alpha/beta diversity profiling, 
             multivariable clinical tracking, and host-pathogen correlation matrices.
"""

import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import matplotlib
matplotlib.use('svg')  # Native vector backend for editable publication-quality figures

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from matplotlib.patches import Patch

# =============================================================================
# GLOBAL REPRODUCIBILITY & STYLE SETUP (NATURE / ISME JOURNAL COMPLIANT)
# =============================================================================
np.random.seed(42)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'  # Retains texts as interactive/editable vectors
sns.set_context("paper", font_scale=1.1)
sns.set_style("white")

print("[1/7] Graphics backend and global style properties successfully initialized.")

# =============================================================================
# COHORT DATA LOADING & UNCOUPLED METADATA HARMONIZATION
# =============================================================================
print("\n[2/7] Loading and harmonizing multi-omic cohorts...")

# 1. Load master metadata and sanitize string variables
df_meta_master = pd.read_csv("metadatos_enriquecidos.csv")
df_meta_master.columns = df_meta_master.columns.str.strip()

for col in df_meta_master.columns:
    if df_meta_master[col].dtype == 'object':
        df_meta_master[col] = df_meta_master[col].astype(str).str.replace('"', '').str.strip()

df_meta_master['Sample_Name'] = df_meta_master['Sample_Name'].astype(str)

# Map clinical covariates generated during dataframe merges (.x suffix resolution)
parasite_cols = {
    'Region.x': 'Region', 
    'Giardia.x': 'Giardia', 
    'Blasto.x': 'Blasto', 
    'EH.x': 'EH', 
    'Crypto.x': 'Crypto', 
    'Ascaris.x': 'Ascaris', 
    'Trichuris.x': 'Trichuris',
    'Polyparasitism.x': 'Polyparasitism'
}
df_meta_master = df_meta_master.rename(columns=parasite_cols)

# Failsafe validation for key categorical variables
if 'Polyparasitism' not in df_meta_master.columns:
    raise KeyError(
        "Critical covariate 'Polyparasitism' or 'Polyparasitism.x' missing from metadata. "
        "Verify column definitions in 'metadatos_enriquecidos.csv'."
    )

df_meta_master['Polyparasitism'] = (
    df_meta_master['Polyparasitism'].astype(str).str.strip().str.capitalize()
)

# 2. MetaPhlAn Cohort Harmonization (Taxonomic Profiling)
df_meta4 = pd.read_csv("tabla_combinada_metaphlan.txt", sep=",", skiprows=1)
df_meta4.columns = df_meta4.columns.str.strip()
df_meta4 = df_meta4.rename(columns={'clade_name': 'Feature'})

samples_metaphlan_raw = set(df_meta4.columns) - {'Feature'}
samples_tax_harmonized = sorted(list(set(df_meta_master['Sample_Name']).intersection(samples_metaphlan_raw)))

# 3. HUMAnN Cohort Harmonization (Functional Profiling)
df_meta_raw = pd.read_csv("metadata.csv")
df_meta_raw.columns = df_meta_raw.columns.str.strip()
df_meta_raw = df_meta_raw.dropna(subset=['Sample_Name', 'humann_id'])
df_meta_raw['Sample_Name'] = df_meta_raw['Sample_Name'].astype(str).str.strip()
df_meta_raw['humann_id'] = df_meta_raw['humann_id'].astype(str).str.strip()

dict_humann_to_sample = pd.Series(df_meta_raw.Sample_Name.values, index=df_meta_raw.humann_id).to_dict()

df_humann = pd.read_csv("humann_merged_all_samples_final_renamed_from_map.tsv", sep="\t")
df_humann = df_humann.rename(columns={'Gene Family': 'Feature'})
df_humann.columns = df_humann.columns.str.strip().str.replace(r'\s+', '', regex=True)

df_humann = df_humann.rename(columns=dict_humann_to_sample)
samples_humann_raw = set(df_humann.columns) - {'Feature'}
samples_func_harmonized = sorted(list(set(df_meta_master['Sample_Name']).intersection(samples_humann_raw)))

print("=" * 60)
print("             COHORT HARMONIZATION INTEGRITY REPORT")
print("=" * 60)
print(f"-> MetaPhlAn Harmonized Cohort (Taxonomy) : {len(samples_tax_harmonized)} samples")
print(f"-> HUMAnN Harmonized Cohort (Functions)   : {len(samples_func_harmonized)} samples")
print("=" * 60)

# =============================================================================
# FIGURE 1: TAXONOMIC STACKED BARPLOT (CORE MICROBIOME VARIABILITY)
# =============================================================================
print("\n[3/7] Generating Figure 1 (Taxonomic Composition)...")

df_species = df_meta4[df_meta4['Feature'].str.contains('s__') & ~df_meta4['Feature'].str.contains('t__')].copy()
df_species['Feature'] = df_species['Feature'].apply(lambda x: x.split('|')[-1].replace('s__', '').replace('_', ' '))
df_species = df_species.set_index('Feature')

df_species_filtered = df_species[samples_tax_harmonized].fillna(0)
df_meta_tax = df_meta_master[df_meta_master['Sample_Name'].isin(samples_tax_harmonized)].set_index('Sample_Name').loc[samples_tax_harmonized]

# Extract and re-normalize the top dominant species to resolve structural composition variations
top_15_taxa = df_species_filtered.mean(axis=1).nlargest(15).index
df_top_taxa = df_species_filtered.loc[top_15_taxa].copy()
df_grouped_taxa = (df_top_taxa / df_top_taxa.sum()) * 100

# Stratify biological samples strictly by geographic origin
df_meta_tax = df_meta_tax.sort_values(by='Region')
df_grouped_taxa = df_grouped_taxa[df_meta_tax.index]

fig, ax = plt.subplots(figsize=(14, 6.5))
df_grouped_taxa.T.plot(kind='bar', stacked=True, cmap='tab20', ax=ax, width=0.85, edgecolor='black', linewidth=0.2)

ax.set_title('Relative Abundance Distribution Within Top 15 Dominant Species (Stratified by Region)', fontsize=13, fontweight='bold', pad=15)
ax.set_ylabel('Within-Top-Taxa Relative Abundance (%)', fontsize=11)
ax.set_xlabel('Samples (Grouped by Geographical Region)', fontsize=11)
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False, fontsize=9.5, prop={'style': 'italic'})

sns.despine()
plt.tight_layout()
plt.savefig('Figure_1_Taxonomic_Abundance_Barplot.svg', format='svg', bbox_inches='tight')

# =============================================================================
# FIGURE 2: HIERARCHICAL CLUSTERMAP WITH CLINICAL COVARIATE TRACKING
# =============================================================================
print("\n[4/7] Generating Figure 2 (Hierarchical Heatmap & Covariates)...")

df_heatmap_data = np.log1p(df_species_filtered.loc[top_15_taxa])
df_anno_tax = df_meta_master[df_meta_master['Sample_Name'].isin(df_heatmap_data.columns)].set_index('Sample_Name').loc[df_heatmap_data.columns]

target_variables = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']

# Define publication-standard discrete palettes for infectious status tracking
color_map_pos_neg = {'POS': '#e41a1c', 'NEG': '#f0f0f0'}
color_map_polyparasitism = {'Poly': '#e41a1c', 'Mono': '#ff7f00', 'Neg': '#f0f0f0'}

unique_regions = df_anno_tax['Region'].dropna().unique()
region_colors = dict(zip(unique_regions, sns.color_palette('Set2', len(unique_regions))))

df_colors_matrix = pd.DataFrame(index=df_anno_tax.index)
for var in target_variables:
    if var == 'Region':
        df_colors_matrix[var] = df_anno_tax[var].map(region_colors)
    elif var == 'Polyparasitism':
        df_colors_matrix[var] = df_anno_tax[var].map(color_map_polyparasitism).fillna('#f0f0f0')
    else:
        df_colors_matrix[var] = df_anno_tax[var].map(color_map_pos_neg).fillna('#f0f0f0')

g = sns.clustermap(
    df_heatmap_data, metric='euclidean', method='average', cmap='viridis',
    col_colors=df_colors_matrix, figsize=(15, 11),
    cbar_kws={'label': 'Log-transformed Abundance log1p(%)'},
    yticklabels=True, xticklabels=False
)

g.ax_heatmap.set_ylabel('Top 15 Dominant Bacterial Species', fontsize=12, fontweight='bold')
g.ax_heatmap.set_xlabel('Study Samples (MetaPhlAn Cohort)', fontsize=12, fontweight='bold')
g.fig.suptitle('Hierarchical Clustering of Taxonomic Profiles with Multivariable Metadata Tracking', fontsize=14, fontweight='bold', y=1.02)

# Build integrated metadata legend map
legend_elements = [Patch(facecolor='black', alpha=0, label='=== Region ===')]
for reg, clr in region_colors.items():
    legend_elements.append(Patch(facecolor=clr, label=f"Region: {reg}"))
legend_elements.append(Patch(facecolor='black', alpha=0, label='=== Parasite Tracking ==='))
legend_elements.append(Patch(facecolor='#e41a1c', label='Positive (POS)'))
legend_elements.append(Patch(facecolor='#f0f0f0', label='Negative (NEG)'))
legend_elements.append(Patch(facecolor='black', alpha=0, label='=== Polyparasitism ==='))
legend_elements.append(Patch(facecolor=color_map_polyparasitism['Poly'], label='Poly (Multiple parasites)'))
legend_elements.append(Patch(facecolor=color_map_polyparasitism['Mono'], label='Mono (Single parasite)'))
legend_elements.append(Patch(facecolor=color_map_polyparasitism['Neg'], label='Neg (No parasites)'))

g.ax_heatmap.legend(handles=legend_elements, bbox_to_anchor=(1.20, 1.01), loc='upper left', frameon=False, fontsize=10)
plt.savefig('Figure_2_Taxonomic_Hierarchical_Heatmap.svg', format='svg', bbox_inches='tight')

# =============================================================================
# FIGURE 3: ALPHA DIVERSITY MULTI-PANEL EVALUATION (3x3 GRID)
# =============================================================================
print("\n[5/7] Generating Figure 3 (Alpha Diversity 3x3 Facets)...")

df_alpha_tax = df_meta_master[df_meta_master['Sample_Name'].isin(samples_tax_harmonized)].dropna(subset=['Shannon'])
analysis_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']

fig, axes = plt.subplots(3, 3, figsize=(18, 15))
axes = axes.flatten()

for i, var in enumerate(analysis_vars):
    ax = axes[i]
    df_sub = df_alpha_tax.dropna(subset=[var])
    groups_list = df_sub[var].unique()
    
    # Select non-parametric test depending on group density
    if len(groups_list) == 2:
        g1 = df_sub[df_sub[var] == groups_list[0]]['Shannon']
        g2 = df_sub[df_sub[var] == groups_list[1]]['Shannon']
        stat, p_val = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        stat_title = f"Mann-Whitney p={p_val:.4f}"
    elif len(groups_list) > 2:
        group_data = [df_sub[df_sub[var] == g]['Shannon'].values for g in groups_list]
        stat, p_val = stats.kruskal(*group_data)
        stat_title = f"Kruskal-Wallis p={p_val:.4f}"
    else:
        stat_title = "Single Group"
        
    sns.boxplot(data=df_sub, x=var, y='Shannon', palette='Set2', width=0.4, ax=ax, fliersize=0)
    sns.stripplot(data=df_sub, x=var, y='Shannon', color='.3', alpha=0.35, size=4, jitter=0.15, ax=ax)
    
    ax.set_title(f"Alpha Diversity by {var}\n({stat_title})", fontsize=11, fontweight='bold', color='navy')
    ax.set_ylabel('Shannon Diversity Index' if i % 3 == 0 else '')
    ax.set_xlabel(f'{var} Status')
    sns.despine(ax=ax)

# Panel 9: Cohort quality check and metadata sample distribution
ax_last = axes[-1]
sns.countplot(data=df_alpha_tax, x='Region', palette='Set2', ax=ax_last, edgecolor='black', linewidth=0.5)
ax_last.set_title("Cohort Sample Distribution\n(Metadata Summary)", fontsize=11, fontweight='bold', color='darkgreen')
ax_last.set_ylabel("Number of Samples")
ax_last.set_xlabel("Geographical Region")

for p in ax_last.patches:
    ax_last.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center', xytext=(0, 5), textcoords='offset points', fontsize=9, fontweight='bold')
sns.despine(ax=ax_last)

plt.tight_layout()
plt.savefig('Figure_3_Alpha_Diversity_Taxonomic_Cohort.svg', format='svg', bbox_inches='tight')

# =============================================================================
# FIGURE 4: BETA DIVERSITY INTEGRATED ORDINATION & NON-PARAMETRIC TESTING
# =============================================================================
print("\n[5b/7] Generating Figure 4 (Beta Diversity Facets with PERMANOVA & PERMDISP)...")

# Principal Component Analysis on log-transformed relative abundances
df_beta_input = np.log1p(df_species_filtered.T)
pca_beta = PCA(n_components=2)
pca_coordinates = pca_beta.fit_transform(df_beta_input)

df_pca_coords = pd.DataFrame(data=pca_coordinates, columns=['PC1', 'PC2'], index=df_beta_input.index).join(df_meta_tax)
var_pc1 = pca_beta.explained_variance_ratio_[0] * 100
var_pc2 = pca_beta.explained_variance_ratio_[1] * 100

beta_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']

fig, axes = plt.subplots(3, 3, figsize=(22, 15), sharex=True, sharey=True)
axes = axes.flatten()

color_pos_neg = {'POS': '#e41a1c', 'NEG': '#b3cde3'}
color_poly = {'Poly': '#e41a1c', 'Mono': '#ff7f00', 'Neg': '#b3cde3'}
color_regions = dict(zip(df_pca_coords['Region'].dropna().unique(), sns.color_palette('Set2', len(df_pca_coords['Region'].dropna().unique()))))

n_perms = 999

for i, var in enumerate(beta_vars):
    ax = axes[i]
    df_sub_plot = df_pca_coords.dropna(subset=[var])

    if var == 'Region':
        current_palette = color_regions
    elif var == 'Polyparasitism':
        current_palette = color_poly
    else:
        current_palette = color_pos_neg
    
    # Real-time mathematical computation of PERMANOVA (Sums of Squares) per panel
    X_stat = df_sub_plot[['PC1', 'PC2']].values
    groups_stat = df_sub_plot[var].values
    unique_groups = np.unique(groups_stat)
    
    if len(unique_groups) >= 2:
        dist_matrix = squareform(pdist(X_stat, metric='euclidean'))
        n_samples = len(groups_stat)
        sst = np.sum(dist_matrix**2) / (2 * n_samples)
        
        ssw = 0.0
        for g in unique_groups:
            g_idx = np.where(groups_stat == g)[0]
            if len(g_idx) > 0:
                ssw += np.sum(dist_matrix[g_idx[:, None], g_idx]**2) / (2 * len(g_idx))
        ssa = sst - ssw
        df_between = len(unique_groups) - 1
        df_within = n_samples - len(unique_groups)
        pseudo_f = (ssa / df_between) / (ssw / df_within) if ssw > 0 else 0
        
        # Empirical permutations for PERMANOVA p-value estimation
        perm_counts = 0
        for _ in range(n_perms):
            p_groups = np.random.permutation(groups_stat)
            p_ssw = 0.0
            for g in unique_groups:
                g_idx = np.where(p_groups == g)[0]
                if len(g_idx) > 0:
                    p_ssw += np.sum(dist_matrix[g_idx[:, None], g_idx]**2) / (2 * len(g_idx))
            p_f = ((sst - p_ssw) / df_between) / (p_ssw / df_within) if p_ssw > 0 else 0
            if p_f >= pseudo_f:
                perm_counts += 1
        permanova_p = (perm_counts + 1) / (n_perms + 1)
        
        # PERMDISP analysis: distances to group centroids followed by parametric ANOVA
        centroids = {g: np.mean(X_stat[groups_stat == g], axis=0) for g in unique_groups}
        dist_to_c = [np.linalg.norm(X_stat[idx] - centroids[g]) for idx, g in enumerate(groups_stat)]
        disp_groups = [np.array(dist_to_c)[groups_stat == g] for g in unique_groups]
        _, permdisp_p = stats.f_oneway(*disp_groups)
        
        stat_label = f"PERMANOVA p={permanova_p:.4f}\nPERMDISP p={permdisp_p:.4f}"
    else:
        stat_label = "Insufficient Groups"

    sns.scatterplot(
        data=df_sub_plot, x='PC1', y='PC2', hue=var, palette=current_palette,
        s=75, alpha=0.85, edgecolor='black', linewidth=0.5, ax=ax
    )
    
    ax.set_title(f'Stratified by {var}', fontsize=12, fontweight='bold', pad=12, color='dimgray')
    
    # Standard transparent bounding text box for statistical readouts
    ax.text(0.05, 0.05, stat_label, transform=ax.transAxes, fontsize=9.5,
            fontweight='semibold', color='darkslategray',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85, edgecolor='none'))
    
    ax.legend(title=var, loc='upper right', frameon=True, framealpha=0.8, edgecolor='none', fontsize=9, title_fontsize=9)
    
    if i >= 6: ax.set_xlabel(f'PC1 ({var_pc1:.1f}%)', fontsize=10.5)
    if i % 3 == 0: ax.set_ylabel(f'PC2 ({var_pc2:.1f}%)', fontsize=10.5)
    sns.despine(ax=ax)

# Panel 9: Global Structural Density Estimation (KDE) to anchor grid symmetry
ax_last = axes[-1]
sns.kdeplot(data=df_pca_coords, x='PC1', y='PC2', cmap='Blues', fill=True, thresh=0.05, ax=ax_last, alpha=0.6)
ax_last.set_title('Global Structural Density\n(All Samples Contours)', fontsize=12, fontweight='bold', pad=10, color='navy')
ax_last.set_xlabel(f'PC1 ({var_pc1:.1f}%)', fontsize=10.5)
sns.despine(ax=ax_last)

fig.suptitle('Beta Diversity Facet Panel - Community Structure (MetaPhlAn Taxonomy Cohort)', fontsize=15, fontweight='bold', y=0.99)
plt.tight_layout()
plt.savefig('Figure_4_Beta_Diversity_Facets.svg', format='svg', bbox_inches='tight')

# =============================================================================
# STATISTICAL EXPORT: STANDALONE PERMANOVA & PERMDISP EVALUATIONS
# =============================================================================
print("\n[5c/7] Exporting formal Beta Diversity statistical matrices...")

def calculate_permanova_and_dispersion(data_matrix, grouping_series, n_permutations=999):
    valid_samples = grouping_series.dropna().index
    X = data_matrix.loc[valid_samples].values
    groups = grouping_series.loc[valid_samples].values
    unique_groups = np.unique(groups)
    
    if len(unique_groups) < 2:
        return "Insufficient groups", "Insufficient groups"
        
    dist_matrix = squareform(pdist(X, metric='euclidean'))
    n_samples = len(groups)
    
    # Sum of Squares Total (SST)
    sst = np.sum(dist_matrix**2) / (2 * n_samples)
    
    # Sum of Squares Within (SSW)
    ssw = 0.0
    for g in unique_groups:
        g_idx = np.where(groups == g)[0]
        if len(g_idx) > 0:
            ssw += np.sum(dist_matrix[g_idx[:, None], g_idx]**2) / (2 * len(g_idx))
            
    ssa = sst - ssw
    df_between = len(unique_groups) - 1
    df_within = n_samples - len(unique_groups)
    pseudo_f_real = (ssa / df_between) / (ssw / df_within) if ssw > 0 else 0
    
    perm_f_counts = 0
    for _ in range(n_permutations):
        perm_groups = np.random.permutation(groups)
        perm_ssw = 0.0
        for g in unique_groups:
            g_idx = np.where(perm_groups == g)[0]
            if len(g_idx) > 0:
                perm_ssw += np.sum(dist_matrix[g_idx[:, None], g_idx]**2) / (2 * len(g_idx))
        perm_ssa = sst - perm_ssw
        perm_f = (perm_ssa / df_between) / (perm_ssw / df_within) if perm_ssw > 0 else 0
        if perm_f >= pseudo_f_real:
            perm_f_counts += 1
            
    permanova_p = (perm_f_counts + 1) / (n_permutations + 1)
    permanova_res = f"Pseudo-F: {pseudo_f_real:.3f} (p={permanova_p:.4f})"
    
    # Beta dispersion centroids tracking
    centroids = {}
    distances_to_centroid = []
    for g in unique_groups:
        g_idx = np.where(groups == g)[0]
        centroids[g] = np.mean(X[g_idx], axis=0)
        
    for i, g in enumerate(groups):
        dist_c = np.linalg.norm(X[i] - centroids[g])
        distances_to_centroid.append(dist_c)
        
    df_disp = pd.DataFrame({'Group': groups, 'DistToCentroid': distances_to_centroid})
    disp_means = df_disp.groupby('Group')['DistToCentroid'].mean().to_dict()
    
    disp_groups_data = [df_disp[df_disp['Group'] == g]['DistToCentroid'].values for g in unique_groups]
    f_stat, permdisp_p = stats.f_oneway(*disp_groups_data)
    
    disp_res = f"F-stat: {f_stat:.3f} (p={permdisp_p:.4f}) -> Means: " + ", ".join([f"{k}={v:.3f}" for k, v in disp_means.items()])
    return permanova_res, disp_res

df_beta_matrix = pd.DataFrame(data=pca_coordinates, index=df_beta_input.index)
stats_summary = []

for var in beta_vars:
    permanova_out, betadisp_out = calculate_permanova_and_dispersion(df_beta_matrix, df_pca_coords[var])
    stats_summary.append({
        'Variable': var,
        'PERMANOVA_Result': permanova_out,
        'Beta_Dispersion_Result': betadisp_out
    })

df_beta_stats = pd.DataFrame(stats_summary)
df_beta_stats.to_csv("Table_2_Beta_Diversity_Statistical_Validation.csv", index=False)
print("[SUCCESS] Statistical verification matrices exported to 'Table_2_Beta_Diversity_Statistical_Validation.csv'.")

# =============================================================================
# UNCOUPLED METABOLIC PATHWAY ANALYSIS (HUMAnN DIFF ABUNDANCE BY REGION)
# =============================================================================
print("\n[6/7] Running functional differential abundance screening (HUMAnN)...")

df_functional = df_humann[~df_humann['Feature'].str.contains(r'\|g__', regex=True)].copy()
df_functional = df_functional.set_index('Feature')

df_func_filtered = df_functional[samples_func_harmonized].fillna(0)
df_meta_func = df_meta_master[df_meta_master['Sample_Name'].isin(samples_func_harmonized)].set_index('Sample_Name').loc[samples_func_harmonized]
df_func_filtered = df_func_filtered[df_func_filtered.std(axis=1) > 0]

region_col_diff = 'Region'
target_regions = df_meta_func[region_col_diff].dropna().unique()

if len(target_regions) >= 2:
    rA_samples = df_meta_func[df_meta_func[region_col_diff] == target_regions[0]].index
    rB_samples = df_meta_func[df_meta_func[region_col_diff] == target_regions[1]].index
    
    p_vals, features, m_A, m_B = [], [], [], []
    
    for feature, row in df_func_filtered.iterrows():
        stat, p = stats.mannwhitneyu(row[rA_samples], row[rB_samples], alternative='two-sided')
        p_vals.append(p)
        features.append(feature)
        m_A.append(row[rA_samples].mean())
        m_B.append(row[rB_samples].mean())
        
    df_diff_res = pd.DataFrame({
        'Functional_Feature': features,
        f'Mean_{target_regions[0]}': m_A,
        f'Mean_{target_regions[1]}': m_B,
        'p_value': p_vals
    }).sort_values('p_value')
    
    # Benjamini-Hochberg False Discovery Rate (FDR) Multi-test Correction
    m = len(df_diff_res)
    df_diff_res['FDR_q_value'] = df_diff_res['p_value'] * m / np.arange(1, m + 1)
    df_diff_res['FDR_q_value'] = np.minimum.accumulate(df_diff_res['FDR_q_value'][::-1])[::-1]
    df_diff_res['FDR_q_value'] = df_diff_res['FDR_q_value'].clip(upper=1.0)
    
    df_sig = df_diff_res[df_diff_res['FDR_q_value'] < 0.05]
    df_sig.to_csv("Table_1_Differential_Abundance_HUMAnN_by_Region.csv", index=False)
    
    print("=" * 60)
    print("         FUNCTIONAL UNCOUPLED MULTI-TEST REGION REPORT")
    print("=" * 60)
    print(f"-> Total Metabolic Gene Families Processed : {m}")
    print(f"-> Significant Features by Region (FDR < 0.05): {len(df_sig)}")
    print("=" * 60)
else:
    print("[ERROR] Insufficient region factors to run functional differential testing.")

# =============================================================================
# FIGURE 5: PATHOGEN LOAD VS HOST ALPHA DIVERSITY CORRELATION DIAGNOSTICS
# =============================================================================
print("\n[7/7] Generating Figure 5 (Pathogen Load vs Diversity Diagnostics Matrix)...")

name_map = {
    'Giardia.load': 'Giardia',
    'Blasto.load': 'Blastocystis',
    'EH.load': 'Entamoeba histolytica',
    'Crypto.load': 'Cryptosporidium',
    'Ascaris.load': 'Ascaris',
    'Trichuris.load': 'Trichuris',
    'Strongy.load': 'Strongyloides'
}
div_name_map = {
    'Shannon': 'Shannon',
    'Observed_Richness': 'Observed Richness',
    'Simpson': 'Simpson',
    'Pielou_Evenness': 'Pielou Evenness'
}
load_cols = list(name_map.keys())
div_cols = list(div_name_map.keys())

missing_cols = [c for c in load_cols + div_cols if c not in df_meta_master.columns]
if missing_cols:
    raise KeyError(f"Missing clinical variables in master metadata registry: {missing_cols}")

# Generate non-parametric Spearman rank correlation matrix and matched p-values
corr_matrix = df_meta_master[load_cols + div_cols].corr(method='spearman')
heatmap_data = corr_matrix.loc[load_cols, div_cols]

pval_matrix = pd.DataFrame(index=load_cols, columns=div_cols, dtype=float)
for lc in load_cols:
    for dc in div_cols:
        _, p_val = stats.spearmanr(df_meta_master[lc], df_meta_master[dc], nan_policy='omit')
        pval_matrix.loc[lc, dc] = p_val

def sig_stars(p):
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return ''

annot_labels = heatmap_data.copy().astype(object)
for lc in load_cols:
    for dc in div_cols:
        annot_labels.loc[lc, dc] = f"{heatmap_data.loc[lc, dc]:.2f}{sig_stars(pval_matrix.loc[lc, dc])}"

heatmap_data.index = [name_map[c] for c in heatmap_data.index]
heatmap_data.columns = [div_name_map[c] for c in heatmap_data.columns]
annot_labels.index = heatmap_data.index
annot_labels.columns = heatmap_data.columns

# Figure Layout: Integrated heatmap + log1p regression metrics + distribution summaries
fig = plt.figure(figsize=(20, 22))
gs = fig.add_gridspec(4, 3, wspace=0.3, hspace=0.4)

# A. Correlation Heatmap
ax_hm = fig.add_subplot(gs[0, :])
sns.heatmap(heatmap_data, annot=annot_labels, fmt='', cmap='RdBu_r', center=0,
            vmin=-1, vmax=1, ax=ax_hm, annot_kws={"size": 10}, linewidths=0.4, linecolor='white')
ax_hm.set_title('Correlation Heatmap: All Pathogens vs All Diversity Metrics\n(* p<0.05, ** p<0.01, *** p<0.001 — Spearman, uncorrected)',
                fontsize=13, fontweight='bold')

# B. Multi-panel Regressions (Pathogen Load vs Shannon Diversity Index)
for i, col in enumerate(load_cols):
    row = (i // 3) + 1
    col_idx = i % 3
    ax = fig.add_subplot(gs[row, col_idx])

    load_transformed = np.log1p(df_meta_master[col])
    rho, pval = stats.spearmanr(load_transformed, df_meta_master['Shannon'], nan_policy='omit')

    sns.regplot(x=load_transformed, y=df_meta_master['Shannon'],
                ax=ax, color='teal', scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})

    ax.set_title(f"{name_map[col]} vs Shannon", fontsize=12, fontweight='bold')
    ax.set_xlabel('Log1p(Load)')
    ax.set_ylabel('Shannon Index')

    ax.text(0.05, 0.95, f"rho: {rho:.2f}\np: {pval:.2e}",
            transform=ax.transAxes, verticalalignment='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    sns.despine(ax=ax)

# C. Cohort Pathogen Prevalence Profile
prevalence = {}
for col in load_cols:
    valid = df_meta_master[col].dropna()
    prevalence[name_map[col]] = (valid > 0).mean() * 100 if len(valid) > 0 else np.nan
prevalence_series = pd.Series(prevalence).sort_values(ascending=False)

ax_prev = fig.add_subplot(gs[3, 1])
sns.barplot(x=prevalence_series.values, y=prevalence_series.index, hue=prevalence_series.index,
            palette='rocket', legend=False, ax=ax_prev, edgecolor='black', linewidth=0.4)
ax_prev.set_title('Pathogen Prevalence\n(% Samples with Detectable Load)', fontsize=12, fontweight='bold', color='darkgreen')
ax_prev.set_xlabel('Prevalence (%)')
ax_prev.set_ylabel('')

for p in ax_prev.patches:
    width = p.get_width()
    ax_prev.annotate(f'{width:.1f}%', (width, p.get_y() + p.get_height() / 2.),
                     va='center', ha='left', fontsize=9, fontweight='bold', xytext=(4, 0), textcoords='offset points')
sns.despine(ax=ax_prev)

# D. Global Pathogen Load Intensity Distribution (Log1p Transformed)
df_melted = np.log1p(df_meta_master[load_cols]).rename(columns=name_map).melt()
df_melted = df_melted[df_melted['value'] > 0]  # Isolate active profiles for visualization clarity

ax_dist = fig.add_subplot(gs[3, 2])
if not df_melted.empty:
    sns.boxplot(data=df_melted, x='value', y='variable', hue='variable', palette='mako', legend=False, ax=ax_dist, fliersize=1)
ax_dist.set_title('Pathogen Load Intensity\n(Log1p Transformed, Non-Zero)', fontsize=12, fontweight='bold', color='purple')
ax_dist.set_xlabel('Log1p(Load)')
ax_dist.set_ylabel('')
sns.despine(ax=ax_dist)

plt.suptitle('Figure 5: Pathogen Load vs. Host Microbiome Alpha Diversity Diagnostics', fontsize=16, fontweight='bold', y=0.99)
plt.savefig('Figure_5_Pathogen_Load_Alpha_Correlation.svg', format='svg', bbox_inches='tight')
print("[SUCCESS] Pipeline completed successfully. Output plots saved to the working directory.")

```
