#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================================
MASTER SCRIPT FOR INTEGRATED METAGENOMIC AND MICROBIOME ANALYSIS
===============================================================================

This script performs a complete analysis of a metagenomic cohort study,
including taxonomic profiling (MetaPhlAn4), functional profiling (HUMAnN3),
alpha/beta diversity, differential abundance, LEfSe‑like biomarker discovery,
and pathogen‑load correlation analysis.

All statistical procedures are rigorously applied with appropriate corrections:
- CLR transformation for compositional beta‑diversity (Aitchison 1986)
- PERMANOVA + PERMDISP with 999 permutations
- FDR correction (Benjamini‑Hochberg) for multiple testing
- Effect‑size filtering (LDA ≥ 1.0) for LEfSe
- Minimum prevalence threshold (≥5 positives) for pathogen correlations

Figures are exported in PNG, PDF, and SVG formats with publication‑grade styling.

Author:  Sergio CAstaneda
Date:    2026‑08‑11
Version: 2.0 
================================================================================
"""

import os
import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
from matplotlib import gridspec
import re

# -----------------------------------------------------------------------------
# 1.  GLOBAL CONFIGURATION – PUBLICATION‑READY STYLING
# -----------------------------------------------------------------------------
# Use non‑interactive backend for server/headless environments
import matplotlib
matplotlib.use('Agg')

# Nature / ISME Journal style guidelines
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'          # Text remains editable
plt.rcParams['pdf.fonttype'] = 42              # TrueType fonts for PDF
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 12

sns.set_context("paper", font_scale=1.1)
sns.set_style("white")
sns.set_palette("Set2")

# Output directories
DIR_FIGURES = "figures"
DIR_TABLES  = "tables"
DIR_LEFSE   = "LEfSe_Results"
DIR_FUNC    = "Functional_Differential_Results"

for d in [DIR_FIGURES, DIR_TABLES, DIR_LEFSE, DIR_FUNC]:
    os.makedirs(d, exist_ok=True)

print("=" * 80)
print("          INTEGRATED METAGENOMIC ANALYSIS PIPELINE")
print("                   Publication‑Ready Script v2.0")
print("=" * 80)
print(f"Output directories: {DIR_FIGURES}, {DIR_TABLES}, {DIR_LEFSE}, {DIR_FUNC}")
print("=" * 80)

# -----------------------------------------------------------------------------
# 2.  HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def clr_transform(data_df, pseudo_count=1e-6):
    """
    Centered Log‑Ratio (CLR) transformation for compositional data.
    
    CLR is the standard approach for multivariate analysis of compositional
    microbiome data (Aitchison 1986, Gloor et al. 2017). It removes the
    unit‑sum constraint and allows use of Euclidean distances in a
    mathematically principled manner.
    
    Parameters
    ----------
    data_df : pd.DataFrame
        Rows = samples, columns = taxa/features. All values must be ≥ 0.
    pseudo_count : float
        Small constant added to avoid log(0) when zeros are present.
    
    Returns
    -------
    pd.DataFrame
        CLR‑transformed data, same shape as input.
    """
    # Add pseudo‑count to avoid undefined log(0)
    X = data_df + pseudo_count
    logX = np.log(X)
    # Row‑wise geometric mean (centering factor)
    geom_mean = logX.mean(axis=1)
    # CLR: subtract the geometric mean from each log‑transformed value
    clr = logX.sub(geom_mean, axis=0)
    return clr


def clean_functional_label(feature_str, max_len=40):
    """
    Convert long HUMAnN feature names into concise, human‑readable labels.
    
    HUMAnN gene family names often contain reaction identifiers, EC numbers,
    and taxonomic suffixes. This function extracts the core enzyme/function
    name for clean visualisation.
    
    Parameters
    ----------
    feature_str : str
        Raw feature name from HUMAnN output.
    max_len : int
        Maximum length of the returned label (truncated with '…' if exceeded).
    
    Returns
    -------
    str
        Cleaned and shortened label.
    """
    patterns = [
        r':\s*\([^)]*\)\s*([^\[\|]+)',   # Extract after ': (expasy) ...'
        r':\s*([^\[\|]+)',                # Extract after ':'
        r'^([^\[\|]+)'                    # Extract before '[' or '|'
    ]
    for pat in patterns:
        m = re.search(pat, feature_str)
        if m:
            label = m.group(1).strip()
            # Capitalise first letter of each word
            label = ' '.join([w.capitalize() for w in label.split()])
            if len(label) > max_len:
                label = label[:max_len-3] + '…'
            return label
    # Fallback: return first 50 characters
    return feature_str[:50]


def permanova_and_betadisp(X, grouping, n_perm=999):
    """
    Compute PERMANOVA (pseudo‑F + p‑value) and PERMDISP (ANOVA on distances
    to centroid) for a given grouping variable.

    PERMANOVA tests for differences in group centroids; PERMDISP tests for
    homogeneity of multivariate dispersions. Both are essential for
    interpreting beta‑diversity patterns (Anderson 2001, 2006).

    Parameters
    ----------
    X : pd.DataFrame
        Samples × features, with samples as index (already transformed,
        e.g., CLR).
    grouping : pd.Series
        Group labels for each sample (index aligned with X.index).
    n_perm : int
        Number of permutations for PERMANOVA p‑value.

    Returns
    -------
    permanova_res : str
        Formatted result: "Pseudo‑F = X.XXX, p = X.XXXX"
    betadisp_res : str
        Formatted result: "F = X.XXX, p = X.XXXX, mean distances: ..."
    """
    # Remove samples with missing group information
    valid = grouping.dropna().index.intersection(X.index)
    X_sub = X.loc[valid].values
    g = grouping.loc[valid].values
    ug = np.unique(g)

    if len(ug) < 2:
        return "Insufficient groups", "Insufficient groups"

    # --- Euclidean distance matrix ---
    D = squareform(pdist(X_sub, metric='euclidean'))
    n = len(g)
    SST = np.sum(D**2) / (2 * n)

    # Within‑group sum of squares (SSW)
    SSW = 0.0
    for grp in ug:
        idx = np.where(g == grp)[0]
        if len(idx) > 0:
            SSW += np.sum(D[idx[:, None], idx]**2) / (2 * len(idx))
    SSA = SST - SSW
    df_b = len(ug) - 1
    df_w = n - len(ug)
    pseudo_F = (SSA / df_b) / (SSW / df_w) if SSW > 0 else 0

    # Permutation test for PERMANOVA
    perm_F = []
    for _ in range(n_perm):
        perm_g = np.random.permutation(g)
        perm_SSW = 0.0
        for grp in ug:
            idx = np.where(perm_g == grp)[0]
            if len(idx) > 0:
                perm_SSW += np.sum(D[idx[:, None], idx]**2) / (2 * len(idx))
        perm_SSA = SST - perm_SSW
        perm_F.append((perm_SSA / df_b) / (perm_SSW / df_w) if perm_SSW > 0 else 0)
    p_val = (np.sum(np.array(perm_F) >= pseudo_F) + 1) / (n_perm + 1)
    permanova_res = f"Pseudo‑F = {pseudo_F:.3f}, p = {p_val:.4f}"

    # --- PERMDISP: distances to group centroids ---
    centroids = {grp: np.mean(X_sub[g == grp], axis=0) for grp in ug}
    dist_centroid = np.array([
        np.linalg.norm(X_sub[i] - centroids[g[i]]) for i in range(n)
    ])
    groups_disp = [dist_centroid[g == grp] for grp in ug]
    F_stat, p_disp = stats.f_oneway(*groups_disp)
    disp_means = {grp: np.mean(dist_centroid[g == grp]) for grp in ug}
    mean_str = ', '.join([f"{k}={v:.3f}" for k, v in disp_means.items()])
    betadisp_res = f"F = {F_stat:.3f}, p = {p_disp:.4f}  [means: {mean_str}]"

    return permanova_res, betadisp_res


# -----------------------------------------------------------------------------
# 3.  DATA LOADING AND HARMONISATION
# -----------------------------------------------------------------------------
print("\n[1/9] Loading and harmonising datasets...")

# 3.1 Master metadata
df_meta = pd.read_csv("metadatos_enriquecidos.csv")
df_meta.columns = df_meta.columns.str.strip()
for col in df_meta.select_dtypes(include='object').columns:
    df_meta[col] = df_meta[col].astype(str).str.replace('"', '').str.strip()
df_meta['Sample_Name'] = df_meta['Sample_Name'].astype(str)

# Standardise column names (remove '.x' suffix from merge)
rename_dict = {
    'Region.x': 'Region',
    'Giardia.x': 'Giardia',
    'Blasto.x': 'Blasto',
    'EH.x': 'EH',
    'Crypto.x': 'Crypto',
    'Ascaris.x': 'Ascaris',
    'Trichuris.x': 'Trichuris',
    'Polyparasitism.x': 'Polyparasitism'
}
df_meta.rename(columns=rename_dict, inplace=True)
if 'Polyparasitism' not in df_meta.columns:
    raise KeyError("Column 'Polyparasitism' not found. Check metadata.")
df_meta['Polyparasitism'] = df_meta['Polyparasitism'].astype(str).str.capitalize()
print(f"  → Master metadata: {len(df_meta)} samples")

# 3.2 MetaPhlAn4 taxonomic profiles
df_meta4 = pd.read_csv("tabla_combinada_metaphlan.txt", sep=",", skiprows=1)
df_meta4.columns = df_meta4.columns.str.strip()
df_meta4.rename(columns={'clade_name': 'Feature'}, inplace=True)
samples_meta_raw = set(df_meta4.columns) - {'Feature'}
samples_tax = sorted(list(set(df_meta['Sample_Name']).intersection(samples_meta_raw)))
print(f"  → MetaPhlAn4 harmonised samples: {len(samples_tax)}")

# 3.3 HUMAnN3 functional profiles
#    IMPORTANT: We use ONLY unstratified rows (no '|' in Feature) to avoid
#    double‑counting taxonomic contributions (including '|unclassified').
df_humann = pd.read_csv("humann_merged_all_samples_final_renamed_from_map.tsv", sep="\t")
df_humann.rename(columns={'Gene Family': 'Feature'}, inplace=True)
df_humann.columns = df_humann.columns.str.strip().str.replace(r'\s+', '', regex=True)

# Keep only unstratified total abundances (no '|' anywhere)
df_humann_unstrat = df_humann[~df_humann['Feature'].str.contains(r'\|', regex=True)].copy()
df_humann_unstrat.set_index('Feature', inplace=True)

# Map humann_id to Sample_Name
map_df = pd.read_csv("metadata.csv")
map_df.columns = map_df.columns.str.strip()
map_df = map_df.dropna(subset=['Sample_Name', 'humann_id'])
map_df['Sample_Name'] = map_df['Sample_Name'].astype(str).str.strip()
map_df['humann_id'] = map_df['humann_id'].astype(str).str.strip()
id_to_sample = pd.Series(map_df.Sample_Name.values, index=map_df.humann_id).to_dict()

df_humann_unstrat.rename(columns=id_to_sample, inplace=True)

samples_func_raw = set(df_humann_unstrat.columns)
samples_func = sorted(list(set(df_meta['Sample_Name']).intersection(samples_func_raw)))
print(f"  → HUMAnN3 harmonised samples: {len(samples_func)}")

# 3.4 Verify normalisation (CPM expected)
func_data = df_humann_unstrat[samples_func].fillna(0)
sample_sums = func_data.sum()
mean_sum = sample_sums.mean()
if mean_sum > 1e5:
    norm_status = "CPM (already normalised)"
elif mean_sum > 1:
    norm_status = "Relative abundance → converting to CPM"
    func_data = func_data.div(func_data.sum(axis=0), axis=1) * 1e6
else:
    norm_status = "Raw counts → converting to CPM"
    func_data = func_data.div(func_data.sum(axis=0), axis=1) * 1e6

# Remove features with zero variance
func_data = func_data[func_data.std(axis=1) > 0]
df_func = func_data
print(f"  → Functional features (non‑zero variance): {df_func.shape[0]}")
print(f"  → Normalisation: {norm_status}")
print(f"  → Database: UniRef90 (v2021_03)")
print(f"  → Stratification: Unstratified (community‑level totals)")

# 3.5 Extract species‑level data (only s__ not t__)
df_species_raw = df_meta4[
    df_meta4['Feature'].str.contains('s__') & ~df_meta4['Feature'].str.contains('t__')
].copy()
df_species_raw['Feature'] = df_species_raw['Feature'].apply(
    lambda x: x.split('|')[-1].replace('s__', '').replace('_', ' ')
)
df_species_raw.set_index('Feature', inplace=True)
df_species = df_species_raw[samples_tax].fillna(0)

# Metadata for taxonomic and functional cohorts
df_meta_tax = df_meta[df_meta['Sample_Name'].isin(samples_tax)].set_index('Sample_Name').loc[samples_tax]
df_meta_func = df_meta[df_meta['Sample_Name'].isin(samples_func)].set_index('Sample_Name').loc[samples_func]

print("\n" + "=" * 60)
print("HARMONISATION SUMMARY")
print("=" * 60)
print(f"Taxonomic cohort  : {len(samples_tax)} samples, {df_species.shape[0]} species")
print(f"Functional cohort : {len(samples_func)} samples, {df_func.shape[0]} gene families")
print("=" * 60)

# -----------------------------------------------------------------------------
# 4.  FIGURE 1 – STACKED BARPLOT OF TOP 15 SPECIES
# -----------------------------------------------------------------------------
print("\n[2/9] Generating Figure 1: Taxonomic stacked barplot...")

# Identify the 15 most abundant species (core microbiome)
top15 = df_species.mean(axis=1).nlargest(15).index
df_top15 = df_species.loc[top15].copy()

# Re‑normalise within top‑15 to 100 % for clear visualisation
df_top15_norm = (df_top15 / df_top15.sum()) * 100

# Order samples by Region
df_meta_tax_sorted = df_meta_tax.sort_values('Region')
df_top15_norm = df_top15_norm[df_meta_tax_sorted.index]

fig, ax = plt.subplots(figsize=(14, 6.5))
df_top15_norm.T.plot(kind='bar', stacked=True, cmap='tab20', ax=ax,
                     width=0.85, edgecolor='black', linewidth=0.2)
ax.set_title('Relative Abundance Distribution Within Top 15 Dominant Species',
             fontsize=13, fontweight='bold', pad=15)
ax.set_ylabel('Within‑Top‑Taxa Relative Abundance (%)', fontsize=11)
ax.set_xlabel('Samples (grouped by Geographical Region)', fontsize=11)
# Italicise species names in legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc='upper left',
          frameon=False, fontsize=9.5, prop={'style': 'italic'})
sns.despine()
plt.tight_layout()

for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(f'{DIR_FIGURES}/Figure_1_Taxonomic_Barplot.{fmt}',
                format=fmt, bbox_inches='tight', dpi=300)
plt.close()
print("  ✓ Figure 1 saved (PNG, PDF, SVG)")

# -----------------------------------------------------------------------------
# 5.  FIGURE 2 – HIERARCHICAL HEATMAP WITH METADATA
# -----------------------------------------------------------------------------
print("\n[3/9] Generating Figure 2: Hierarchical heatmap...")

df_heat = np.log1p(df_species.loc[top15])   # log1p for visualisation
df_anno = df_meta_tax

# Colour mappings
region_pal = dict(zip(df_anno['Region'].unique(),
                      sns.color_palette('Set2', len(df_anno['Region'].unique()))))
pos_neg = {'POS': '#e41a1c', 'NEG': '#f0f0f0'}
poly_pal = {'Poly': '#e41a1c', 'Mono': '#ff7f00', 'Neg': '#f0f0f0'}

# Build annotation colour matrix
cols_anno = pd.DataFrame(index=df_anno.index)
for var in ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']:
    if var == 'Region':
        cols_anno[var] = df_anno[var].map(region_pal)
    elif var == 'Polyparasitism':
        cols_anno[var] = df_anno[var].map(poly_pal).fillna('#f0f0f0')
    else:
        cols_anno[var] = df_anno[var].map(pos_neg).fillna('#f0f0f0')

g = sns.clustermap(df_heat, metric='euclidean', method='average', cmap='viridis',
                   col_colors=cols_anno, figsize=(15, 11),
                   cbar_kws={'label': 'Log1p abundance (%)'},
                   yticklabels=True, xticklabels=False)
g.ax_heatmap.set_ylabel('Top 15 Dominant Species', fontsize=12, fontweight='bold')
g.ax_heatmap.set_xlabel('Samples', fontsize=12, fontweight='bold')
g.fig.suptitle('Hierarchical Clustering of Taxonomic Profiles with Metadata Tracking',
               fontsize=14, fontweight='bold', y=1.02)

# Legend
leg_elem = [Patch(facecolor='black', alpha=0, label='=== Region ===')]
for reg, clr in region_pal.items():
    leg_elem.append(Patch(facecolor=clr, label=f"Region: {reg}"))
leg_elem.append(Patch(facecolor='black', alpha=0, label='=== Parasite ==='))
leg_elem.append(Patch(facecolor='#e41a1c', label='Positive'))
leg_elem.append(Patch(facecolor='#f0f0f0', label='Negative'))
leg_elem.append(Patch(facecolor='black', alpha=0, label='=== Polyparasitism ==='))
leg_elem.append(Patch(facecolor=poly_pal['Poly'], label='Poly'))
leg_elem.append(Patch(facecolor=poly_pal['Mono'], label='Mono'))
leg_elem.append(Patch(facecolor=poly_pal['Neg'], label='Neg'))
g.ax_heatmap.legend(handles=leg_elem, bbox_to_anchor=(1.20, 1.01),
                    loc='upper left', frameon=False, fontsize=10)

for fmt in ['png', 'pdf', 'svg']:
    g.savefig(f'{DIR_FIGURES}/Figure_2_Heatmap.{fmt}',
              format=fmt, bbox_inches='tight', dpi=300)
plt.close()
print("  ✓ Figure 2 saved (PNG, PDF, SVG)")

# -----------------------------------------------------------------------------
# 6.  FIGURE 3 – ALPHA DIVERSITY PANEL (Shannon)
# -----------------------------------------------------------------------------
print("\n[4/9] Generating Figure 3: Alpha diversity panel...")

alpha_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']
df_alpha = df_meta_tax.dropna(subset=['Shannon'])

fig, axes = plt.subplots(3, 3, figsize=(18, 15))
axes = axes.flatten()

for i, var in enumerate(alpha_vars):
    ax = axes[i]
    sub = df_alpha.dropna(subset=[var])
    groups = sub[var].unique()
    # Choose statistical test based on number of groups
    if len(groups) == 2:
        g1 = sub[sub[var]==groups[0]]['Shannon']
        g2 = sub[sub[var]==groups[1]]['Shannon']
        _, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        label = f"Mann‑Whitney p={p:.4f}"
    elif len(groups) > 2:
        data = [sub[sub[var]==g]['Shannon'].values for g in groups]
        _, p = stats.kruskal(*data)
        label = f"Kruskal‑Wallis p={p:.4f}"
    else:
        label = "Single group"
    sns.boxplot(data=sub, x=var, y='Shannon', palette='Set2',
                width=0.4, ax=ax, fliersize=0)
    sns.stripplot(data=sub, x=var, y='Shannon', color='.3',
                  alpha=0.35, size=4, jitter=0.15, ax=ax)
    ax.set_title(f"{var}\n({label})", fontsize=11, fontweight='bold', color='navy')
    ax.set_ylabel('Shannon Index' if i % 3 == 0 else '')
    ax.set_xlabel('')
    sns.despine(ax=ax)

# Last panel: sample distribution
ax_last = axes[-1]
sns.countplot(data=df_alpha, x='Region', palette='Set2', ax=ax_last, edgecolor='black')
ax_last.set_title("Cohort distribution", fontsize=11, fontweight='bold', color='darkgreen')
ax_last.set_ylabel("Number of samples")
ax_last.set_xlabel("Region")
for p in ax_last.patches:
    ax_last.annotate(f'{int(p.get_height())}',
                     (p.get_x()+p.get_width()/2., p.get_height()),
                     ha='center', va='center', xytext=(0,5),
                     textcoords='offset points')
sns.despine(ax=ax_last)

plt.tight_layout()
for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(f'{DIR_FIGURES}/Figure_3_Alpha_Diversity.{fmt}',
                format=fmt, bbox_inches='tight', dpi=300)
plt.close()
print("  ✓ Figure 3 saved (PNG, PDF, SVG)")

# -----------------------------------------------------------------------------
# 7.  FIGURE 4 – BETA DIVERSITY WITH CLR + PERMANOVA/PERMDISP
# -----------------------------------------------------------------------------
print("\n[5/9] Generating Figure 4: Beta diversity with CLR transformation...")

# CLR transform species abundances (samples × species)
df_clr = clr_transform(df_species.T)

# PCA on CLR data
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(df_clr)
df_pca = pd.DataFrame(data=pca_coords, columns=['PC1','PC2'], index=df_clr.index)
df_pca = df_pca.join(df_meta_tax)
var1, var2 = pca.explained_variance_ratio_[0]*100, pca.explained_variance_ratio_[1]*100

# Variables for beta panel
beta_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']
color_pos_neg = {'POS': '#e41a1c', 'NEG': '#b3cde3'}
color_poly = {'Poly': '#e41a1c', 'Mono': '#ff7f00', 'Neg': '#b3cde3'}
color_regions = dict(zip(df_pca['Region'].dropna().unique(),
                         sns.color_palette('Set2', len(df_pca['Region'].dropna().unique()))))

fig, axes = plt.subplots(3, 3, figsize=(22, 15), sharex=True, sharey=True)
axes = axes.flatten()
beta_stats = []

for i, var in enumerate(beta_vars):
    ax = axes[i]
    sub = df_pca.dropna(subset=[var])
    pal = (color_regions if var == 'Region' else
           color_poly if var == 'Polyparasitism' else color_pos_neg)

    # Compute PERMANOVA and PERMDISP on CLR data
    perm_res, disp_res = permanova_and_betadisp(df_clr, df_meta_tax[var])
    beta_stats.append({'Variable': var, 'PERMANOVA': perm_res, 'PERMDISP': disp_res})

    sns.scatterplot(data=sub, x='PC1', y='PC2', hue=var, palette=pal,
                    s=75, alpha=0.85, edgecolor='black', linewidth=0.5, ax=ax)
    ax.set_title(f'Stratified by {var}', fontsize=12, fontweight='bold', color='dimgray')
    ax.text(0.05, 0.05, f"PERMANOVA: {perm_res}\nPERMDISP: {disp_res}",
            transform=ax.transAxes, fontsize=8.5, fontweight='semibold',
            color='darkslategray',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      alpha=0.85, edgecolor='none'))
    ax.legend(title=var, loc='upper right', frameon=True,
              framealpha=0.8, edgecolor='none', fontsize=8, title_fontsize=8)
    if i >= 6: ax.set_xlabel(f'PC1 ({var1:.1f}%)', fontsize=10.5)
    if i % 3 == 0: ax.set_ylabel(f'PC2 ({var2:.1f}%)', fontsize=10.5)
    sns.despine(ax=ax)

# Last panel: global KDE
ax_last = axes[-1]
sns.kdeplot(data=df_pca, x='PC1', y='PC2', cmap='Blues',
            fill=True, thresh=0.05, ax=ax_last, alpha=0.6)
ax_last.set_title('Global Structural Density (All Samples)', fontsize=12,
                  fontweight='bold', color='navy')
ax_last.set_xlabel(f'PC1 ({var1:.1f}%)', fontsize=10.5)
sns.despine(ax=ax_last)

fig.suptitle('Beta Diversity Facet Panel – CLR‑transformed Compositional Data',
             fontsize=15, fontweight='bold', y=0.99)
plt.tight_layout()

for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(f'{DIR_FIGURES}/Figure_4_Beta_Diversity.{fmt}',
                format=fmt, bbox_inches='tight', dpi=300)
plt.close()

# Save beta statistics table
pd.DataFrame(beta_stats).to_csv(f'{DIR_TABLES}/Table_Beta_Diversity_Statistics.csv', index=False)
print("  ✓ Figure 4 saved (PNG, PDF, SVG)")
print("  ✓ Beta diversity statistics saved to table")

# -----------------------------------------------------------------------------
# 8.  LEfSe‑LIKE ANALYSIS WITH EFFECT SIZE FILTERING (LDA ≥ 1.0)
# -----------------------------------------------------------------------------
print("\n[6/9] Running LEfSe‑like differential abundance analysis...")

def lefse_analysis(var, df_features, df_metadata, out_dir, top_n=10, min_lda=1.0):
    """
    LEfSe‑style biomarker discovery with statistical filtering and effect‑size
    ranking. Only features with FDR < 0.05 and |LDA| ≥ min_lda are retained.
    """
    valid = df_metadata[df_metadata[var].notna()].index.intersection(df_features.columns)
    if len(valid) < 10:
        return None
    y = df_metadata.loc[valid, var]
    groups = y.unique()

    # 1. Univariate filtering (Mann‑Whitney or Kruskal‑Wallis)
    pvals = {}
    for feat in df_features.index:
        data = df_features.loc[feat, valid]
        if np.var(data) < 1e-10:
            continue
        try:
            if len(groups) == 2:
                g1, g2 = data[y==groups[0]], data[y==groups[1]]
                if len(g1) < 2 or len(g2) < 2:
                    continue
                _, p = stats.ranksums(g1, g2)
            else:
                grp_data = [data[y==g] for g in groups]
                if any(len(g) < 2 for g in grp_data):
                    continue
                _, p = stats.kruskal(*grp_data)
            pvals[feat] = p
        except Exception:
            continue
    if len(pvals) < 2:
        return None

    df_p = pd.DataFrame.from_dict(pvals, orient='index', columns=['P_Value'])
    _, fdr, _, _ = multipletests(df_p['P_Value'], method='fdr_bh')
    df_p['FDR'] = fdr
    sig_feats = df_p[df_p['FDR'] < 0.05].index.tolist()
    if len(sig_feats) < 2:
        return None

    # 2. LDA on significant features only
    X = df_features.loc[sig_feats, valid].T
    lda = LDA(solver='lsqr', shrinkage='auto')
    lda.fit(X, y)

    # 3. Extract LDA scores (effect size)
    coef = lda.coef_[0] if len(lda.classes_) == 2 else lda.coef_[0]
    df_lda = pd.DataFrame({
        'Feature': X.columns,
        'LDA_Score': coef,
        'Abs_LDA': np.abs(coef)
    }).merge(df_p, left_on='Feature', right_index=True)

    # 4. Apply effect‑size filter (LDA ≥ min_lda)
    df_lda = df_lda[df_lda['Abs_LDA'] >= min_lda].sort_values('Abs_LDA', ascending=False)
    if df_lda.empty:
        return None

    # 5. Select top_n positive and top_n negative
    top_pos = df_lda[df_lda['LDA_Score'] > 0].head(top_n)
    top_neg = df_lda[df_lda['LDA_Score'] < 0].head(top_n)
    selected = pd.concat([top_pos, top_neg])

    # 6. Save outputs
    df_lda.to_csv(f'{out_dir}/LEfSe_LDA_Scores_{var}.csv', index=False)
    selected.to_csv(f'{out_dir}/LEfSe_Top{top_n}_{var}.csv', index=False)

    # 7. Plot
    plot_data = pd.concat([df_lda.head(top_n), df_lda.tail(top_n)])
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=plot_data, x='LDA_Score', y='Feature',
                palette='RdBu_r', edgecolor='black', ax=ax)
    ax.set_title(f"Discriminative Features: {var}\n(FDR < 0.05, |LDA| ≥ {min_lda})",
                 fontsize=12, fontweight='bold')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
    ax.set_xlabel("LDA Score (Effect Size)")
    ax.set_ylabel("")
    plt.tight_layout()
    for fmt in ['png', 'pdf', 'svg']:
        plt.savefig(f'{out_dir}/LEfSe_Plot_{var}.{fmt}',
                    format=fmt, bbox_inches='tight', dpi=300)
    plt.close()

    return {
        'Variable': var,
        'Total_significant': len(sig_feats),
        'Passed_effect_filter': len(df_lda),
        'Selected_top': len(selected),
        'Groups': groups.tolist()
    }

# Run for all variables
lefse_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']
lefse_summary = []
for var in lefse_vars:
    print(f"  → Processing {var}...")
    res = lefse_analysis(var, df_species, df_meta_tax, DIR_LEFSE, top_n=10, min_lda=1.0)
    if res:
        lefse_summary.append(res)

pd.DataFrame(lefse_summary).to_csv(f'{DIR_TABLES}/Table_LEfSe_Summary.csv', index=False)

# Multi‑panel summary figure (Figure 5)
full_names = {
    'Region': 'Geographic Region',
    'Giardia': 'Giardia lamblia',
    'Blasto': 'Blastocystis',
    'EH': 'Entamoeba histolytica',
    'Crypto': 'Cryptosporidium sp.',
    'Ascaris': 'Ascaris lumbricoides',
    'Trichuris': 'Trichuris trichiura',
    'Polyparasitism': 'Polyparasitism'
}
vars_with_plots = [v for v in lefse_vars if os.path.exists(f'{DIR_LEFSE}/LEfSe_LDA_Scores_{v}.csv')]

if vars_with_plots:
    n_vars = len(vars_with_plots)
    n_cols = 2
    n_rows = (n_vars + 1) // n_cols
    fig = plt.figure(figsize=(16, 5 * n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols, hspace=0.4, wspace=0.3)

    for idx, var in enumerate(vars_with_plots):
        df = pd.read_csv(f'{DIR_LEFSE}/LEfSe_LDA_Scores_{var}.csv')
        df_top = pd.concat([df.head(10), df.tail(10)])
        ax = fig.add_subplot(gs[idx])
        sns.barplot(data=df_top, x='LDA_Score', y='Feature',
                    palette='RdBu_r', edgecolor='black', ax=ax)
        ax.set_title(f"{full_names.get(var, var)}\n(FDR < 0.05, |LDA| ≥ 1.0)",
                     fontsize=13, fontweight='bold')
        ax.axvline(x=0, color='black', linewidth=0.8, linestyle='--')
        ax.set_xlabel("LDA Score")
        ax.set_ylabel("")
        sns.despine(ax=ax)

    plt.tight_layout()
    for fmt in ['png', 'pdf', 'svg']:
        plt.savefig(f'{DIR_FIGURES}/Figure_5_LEfSe_Multipanel.{fmt}',
                    format=fmt, bbox_inches='tight', dpi=300)
    plt.close()
    print("  ✓ Figure 5 saved (PNG, PDF, SVG)")

print("  ✓ LEfSe analysis complete")

# -----------------------------------------------------------------------------
# 9.  FUNCTIONAL DIFFERENTIAL ANALYSIS (HUMAnN unstratified)
# -----------------------------------------------------------------------------
print("\n[7/9] Running functional differential analysis...")

def func_diff_test(df_features, grouping):
    valid = grouping.dropna().index.intersection(df_features.columns)
    groups = grouping.loc[valid]
    ug = np.unique(groups.values)
    if len(ug) < 2:
        return None, ug
    df_sub = df_features[valid]
    pvals, feats = [], []
    means = {g: [] for g in ug}
    for feat, row in df_sub.iterrows():
        data_by_group = [row[groups[groups==g].index] for g in ug]
        if len(ug) == 2:
            _, p = stats.mannwhitneyu(data_by_group[0], data_by_group[1],
                                      alternative='two-sided')
        else:
            _, p = stats.kruskal(*data_by_group)
        pvals.append(p)
        feats.append(feat)
        for g, vals in zip(ug, data_by_group):
            means[g].append(vals.mean())
    res = pd.DataFrame({'Functional_Feature': feats, 'p_value': pvals})
    for g in ug:
        res[f'Mean_{g}'] = means[g]
    res = res.sort_values('p_value').reset_index(drop=True)
    # FDR
    m = len(res)
    res['FDR_q_value'] = res['p_value'] * m / np.arange(1, m+1)
    res['FDR_q_value'] = np.minimum.accumulate(res['FDR_q_value'][::-1])[::-1]
    res['FDR_q_value'] = res['FDR_q_value'].clip(upper=1.0)
    return res, ug

func_vars = ['Region', 'Giardia', 'Blasto', 'EH', 'Crypto', 'Ascaris', 'Trichuris', 'Polyparasitism']
func_summary = []

for var in func_vars:
    print(f"  → Processing {var}...")
    res, ug = func_diff_test(df_func, df_meta_func[var])
    if res is None:
        continue
    sig = res[res['FDR_q_value'] < 0.05]
    test = "Mann‑Whitney U" if len(ug) == 2 else "Kruskal‑Wallis"
    out = f'{DIR_FUNC}/Differential_Abundance_HUMAnN_by_{var}.csv'
    sig.to_csv(out, index=False)
    func_summary.append({
        'Variable': var,
        'Test': test,
        'Groups': ', '.join(ug),
        'Total': len(res),
        'Significant_FDR<0.05': len(sig)
    })
    # Generate boxplot for top 15 features
    if not sig.empty:
        top_feats = sig.head(15)['Functional_Feature'].tolist()
        clean_map = {f: clean_functional_label(f) for f in top_feats}
        data_plot = df_func.loc[top_feats].T
        data_plot.rename(columns=clean_map, inplace=True)
        data_plot = data_plot.merge(df_meta_func[[var]], left_index=True, right_index=True)
        melted = data_plot.melt(id_vars=[var], var_name='Feature', value_name='Abundance')
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.boxplot(data=melted, x='Abundance', y='Feature', hue=var,
                    palette='Set2', width=0.6, ax=ax)
        ax.set_title(f"Top Functional Features: {full_names.get(var, var)}", fontsize=14, fontweight='bold')
        ax.set_xlabel("Relative Abundance (CPM)")
        ax.set_ylabel("")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        for fmt in ['png', 'pdf', 'svg']:
            plt.savefig(f'{DIR_FIGURES}/Figure_6_Functional_{var}.{fmt}',
                        format=fmt, bbox_inches='tight', dpi=300)
        plt.close()

pd.DataFrame(func_summary).to_csv(f'{DIR_TABLES}/Table_Functional_Summary.csv', index=False)
print("  ✓ Functional analysis complete (Figure 6 panels saved)")

# -----------------------------------------------------------------------------
# 10. FIGURE 7 – PATHOGEN LOAD VS DIVERSITY CORRELATIONS
# -----------------------------------------------------------------------------
print("\n[8/9] Generating Figure 7: Pathogen load–diversity correlations...")

# 10.1 Define columns
load_cols = ['Giardia.load', 'Blasto.load', 'EH.load', 'Crypto.load',
             'Ascaris.load', 'Trichuris.load', 'Strongy.load']
div_cols = ['Shannon', 'Observed_Richness', 'Simpson', 'Pielou_Evenness']
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
# Keep only columns present
load_cols = [c for c in load_cols if c in df_meta.columns]
div_cols = [c for c in div_cols if c in df_meta.columns]

# 10.2 Correlation heatmap
corr = df_meta[load_cols + div_cols].corr(method='spearman')
heat = corr.loc[load_cols, div_cols]

pvals = pd.DataFrame(index=load_cols, columns=div_cols)
for lc in load_cols:
    for dc in div_cols:
        _, p = stats.spearmanr(df_meta[lc], df_meta[dc], nan_policy='omit')
        pvals.loc[lc, dc] = p

def stars(p):
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return ''

annot = heat.copy().astype(object)
for lc in load_cols:
    for dc in div_cols:
        annot.loc[lc, dc] = f"{heat.loc[lc, dc]:.2f}{stars(pvals.loc[lc, dc])}"

heat.index = [name_map[c] for c in heat.index]
heat.columns = [div_name_map[c] for c in heat.columns]
annot.index = heat.index
annot.columns = heat.columns

# 10.3 Filter pathogens with ≥5 positive samples for scatter plots
valid_pathogens = []
for lc in load_cols:
    valid = df_meta[lc].dropna()
    n_pos = (valid > 0).sum()
    if n_pos >= 5:
        valid_pathogens.append(lc)
        print(f"  {name_map[lc]}: {n_pos} positive samples (kept)")
    else:
        print(f"  {name_map[lc]}: {n_pos} positive samples (excluded, <5)")

# 10.4 Build composite figure
fig = plt.figure(figsize=(20, 22))
gs = fig.add_gridspec(4, 3, wspace=0.3, hspace=0.4)

# A. Heatmap
ax_hm = fig.add_subplot(gs[0, :])
sns.heatmap(heat, annot=annot, fmt='', cmap='RdBu_r', center=0,
            vmin=-1, vmax=1, ax=ax_hm, annot_kws={"size": 10},
            linewidths=0.4, linecolor='white')
ax_hm.set_title('Spearman Correlations: Pathogen Load vs Diversity Metrics\n(* p<0.05, ** p<0.01, *** p<0.001, uncorrected)',
                fontsize=13, fontweight='bold')

# B. Scatter plots (only for valid pathogens, max 7)
for i, lc in enumerate(valid_pathogens[:7]):
    row = (i // 3) + 1
    col = i % 3
    ax = fig.add_subplot(gs[row, col])
    x = np.log1p(df_meta[lc])
    y = df_meta['Shannon']
    rho, p = stats.spearmanr(x, y, nan_policy='omit')
    sns.regplot(x=x, y=y, ax=ax, color='teal',
                scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
    n_valid = x.notna().sum()
    ax.set_title(f"{name_map[lc]} vs Shannon\n(n={n_valid}, ρ={rho:.2f})",
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('Log1p(Load)')
    ax.set_ylabel('Shannon Index')
    ax.text(0.05, 0.95, f"ρ = {rho:.2f}\np = {p:.2e}",
            transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    sns.despine(ax=ax)

# C. Prevalence barplot with raw numbers
prevalence = {}
for lc in load_cols:
    valid = df_meta[lc].dropna()
    if len(valid) > 0:
        n_pos = (valid > 0).sum()
        prevalence[name_map[lc]] = {
            'percent': (n_pos / len(valid)) * 100,
            'n_pos': n_pos,
            'n_total': len(valid)
        }
prev_df = pd.DataFrame(prevalence).T.sort_values('percent', ascending=False)

ax_prev = fig.add_subplot(gs[3, 1])
bars = ax_prev.barh(prev_df.index, prev_df['percent'],
                    color=sns.color_palette('rocket', len(prev_df)),
                    edgecolor='black', linewidth=0.5)
ax_prev.set_title('Pathogen Prevalence (% Samples with Detectable Load)',
                  fontsize=12, fontweight='bold', color='darkgreen')
ax_prev.set_xlabel('Prevalence (%)')
ax_prev.set_ylabel('')
for i, (idx, row) in enumerate(prev_df.iterrows()):
    ax_prev.text(row['percent'] + 0.5, i,
                 f"{row['n_pos']}/{row['n_total']} ({row['percent']:.1f}%)",
                 va='center', fontsize=9, fontweight='bold')
sns.despine(ax=ax_prev)

# D. Load distribution
melted = np.log1p(df_meta[load_cols]).rename(columns=name_map)
melted = melted.melt(var_name='Pathogen', value_name='Log1p_Load').dropna()
ax_dist = fig.add_subplot(gs[3, 2])
sns.boxplot(data=melted, x='Log1p_Load', y='Pathogen', hue='Pathogen',
            palette='crest', legend=False, ax=ax_dist, width=0.6, fliersize=2)
ax_dist.set_title('Load Distribution Overview (Log1p‑transformed)',
                  fontsize=12, fontweight='bold', color='navy')
ax_dist.set_xlabel('Log1p(Load)')
ax_dist.set_ylabel('')
sns.despine(ax=ax_dist)

fig.suptitle('Pathogen Load – Diversity Correlation Analysis', fontsize=18, fontweight='bold', y=0.985)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

for fmt in ['png', 'pdf', 'svg']:
    plt.savefig(f'{DIR_FIGURES}/Figure_7_Pathogen_Diversity_Correlation.{fmt}',
                format=fmt, bbox_inches='tight', dpi=300)
plt.close()
print("  ✓ Figure 7 saved (PNG, PDF, SVG)")

# -----------------------------------------------------------------------------
# 11. FINAL SUMMARY
# -----------------------------------------------------------------------------
print("\n" + "=" * 80)
print("                              ANALYSIS COMPLETE")
print("=" * 80)
print("\nOutput summary:")
print(f"  Figures  : {DIR_FIGURES}/ (PNG, PDF, SVG)")
print(f"  Tables   : {DIR_TABLES}/ (CSV)")
print(f"  LEfSe    : {DIR_LEFSE}/")
print(f"  Func. diff: {DIR_FUNC}/")
print("\nFigures generated:")
print("  Figure 1: Taxonomic Stacked Barplot (Top 15 Species)")
print("  Figure 2: Hierarchical Heatmap with Metadata Tracking")
print("  Figure 3: Alpha Diversity Panel (Shannon)")
print("  Figure 4: Beta Diversity with CLR + PERMANOVA/PERMDISP")
print("  Figure 5: LEfSe Multi‑panel (Effect‑size filtered)")
print("  Figure 6: Functional Differential Analysis Panels")
print("  Figure 7: Pathogen Load–Diversity Correlations")
print("\nAll statistical procedures are documented within the script.")
print("=" * 80)
