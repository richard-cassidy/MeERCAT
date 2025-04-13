# mercat_analyzer/mercat_analyzer/visualize.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

from src.config import PLOT_HEATMAP_TOP_N_PAIRS
try:
    from adjustText import adjust_text
    adjust_text_available = True
except ImportError:
    adjust_text_available = False
from matplotlib.lines import Line2D

# Import default config values
from .config import PLOT_PADJ_THRESHOLD, PLOT_RHO_THRESHOLD, \
                   PLOT_N_TOP_VOLCANO_LABELS, PLOT_HEATMAP_MAX_FEATURES, \
                   PLOT_N_TOP_POS_NEG_SCATTER

def plot_rho_distribution(df_results, save_path=None):
    """Generates and saves a histogram of Spearman Rho values."""
    print("\n--- Generating Rho Distribution Histogram ---")
    if df_results is None or df_results.empty or 'rho' not in df_results.columns:
        print("Cannot plot: Input data invalid or missing 'rho' column.")
        return

    try:
        plt.figure(figsize=(8, 5))
        sns.histplot(df_results['rho'], bins=50, kde=True)
        plt.title('Distribution of Spearman Correlation Coefficients (Rho)')
        plt.xlabel('Spearman Rho'); plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.5)
        if save_path:
            filename = os.path.join(save_path, "rho_distribution_histogram.png")
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved Rho histogram to: {filename}")
        plt.show(); plt.close()
    except Exception as e: print(f"Error generating Rho histogram: {e}"); plt.close()


def plot_volcano(df_results, p_adj_col='p_adjusted', rho_col='rho',
                 gene_col='gene', metab_col='metabolite',
                 p_adj_thresh=PLOT_PADJ_THRESHOLD, rho_thresh=PLOT_RHO_THRESHOLD,
                 n_labels=PLOT_N_TOP_VOLCANO_LABELS, save_path=None):
    """Generates and saves a labeled Volcano Plot."""
    print("\n--- Generating Volcano Plot ---")
    required_cols = [p_adj_col, rho_col, gene_col, metab_col]
    if df_results is None or df_results.empty or not all(c in df_results.columns for c in required_cols):
        print(f"Cannot plot: Input data invalid or missing required columns: {required_cols}")
        return

    try:
        df_volcano = df_results.copy()
        df_volcano['-log10p'] = -np.log10(df_volcano[p_adj_col].replace(0, 1e-300).fillna(1.0))
        df_volcano['significant'] = (df_volcano[p_adj_col] < p_adj_thresh) & (abs(df_volcano[rho_col]) >= rho_thresh)

        plt.figure(figsize=(13, 9)); ax = plt.gca()
        sns.scatterplot(data=df_volcano, x=rho_col, y='-log10p', hue='significant',
                        palette={True: 'red', False: 'grey'}, alpha=0.6, s=15, legend=False, ax=ax)
        plt.axhline(-np.log10(p_adj_thresh), color='blue', linestyle='--', lw=1)
        plt.axvline(rho_thresh, color='blue', linestyle='--', lw=1); plt.axvline(-rho_thresh, color='blue', linestyle='--', lw=1)

        texts = []
        significant_subset = df_volcano[df_volcano['significant']]
        if not significant_subset.empty:
            top_pos = significant_subset.nlargest(n_labels, rho_col)
            top_neg = significant_subset.nsmallest(n_labels, rho_col)
            points_to_label = pd.concat([top_pos, top_neg]).drop_duplicates()
            print(f"Preparing labels for up to {len(points_to_label)} points...")
            for i, row in points_to_label.iterrows():
                 label_text = f"{row[gene_col]} / {row[metab_col]}"
                 texts.append(plt.text(row[rho_col], row['-log10p'], label_text, fontsize=6))
            if adjust_text_available and texts:
                 print("Adjusting labels...")
                 adjust_text(texts, lim=200, force_points=(0.2, 0.4), force_text=(0.3, 0.6),
                             arrowprops=dict(arrowstyle="-", color='grey', lw=0.5, alpha=0.8))
            elif texts: print("adjustText not available, labels might overlap.")
        else: print("No significant points found for labeling.")

        legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'Sig. (p<{p_adj_thresh}, |rho|>{rho_thresh})', markerfacecolor='red', markersize=6), Line2D([0], [0], marker='o', color='w', label='Not Sig.', markerfacecolor='grey', markersize=6), Line2D([0], [0], color='blue', lw=1, linestyle='--', label='Thresholds')]
        ax.legend(handles=legend_elements, fontsize='medium', loc='upper left', bbox_to_anchor=(1.02, 1))
        plt.title('Volcano Plot: Correlation Strength vs. Significance', fontsize=14)
        plt.xlabel(f'Spearman {rho_col}', fontsize=12); plt.ylabel(f'-log10 ({p_adj_col})', fontsize=12)
        plt.grid(alpha=0.4); plt.tight_layout(rect=[0, 0, 0.82, 1])

        if save_path:
            filename = os.path.join(save_path, "correlation_volcano_plot_labeled.png")
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved Volcano plot to: {filename}")
        plt.show(); plt.close()

    except Exception as e: print(f"Error generating Volcano plot: {e}"); plt.close()


def plot_correlation_clustermap(df_results, rho_col='rho', gene_col='gene', metab_col='metabolite',
                                save_path=None):
    """Generates and saves a clustermap of all correlations."""
    print("\n--- Generating Global Correlation Clustermap ---")
    required_cols = [gene_col, metab_col, rho_col]
    if df_results is None or df_results.empty or not all(c in df_results.columns for c in required_cols):
        print(f"Cannot plot: Input data invalid or missing columns: {required_cols}")
        return

    try:
        print("Pivoting data (may take time)...")
        clustermap_matrix = df_results.pivot_table(index=gene_col, columns=metab_col, values=rho_col)
        nan_count = clustermap_matrix.isna().sum().sum()
        if nan_count > 0: clustermap_matrix.fillna(0, inplace=True); print(f"Filled {nan_count} NaNs with 0.")
        print(f"Matrix shape: {clustermap_matrix.shape}")

        if clustermap_matrix.empty: print("Skipping clustermap: Pivot table is empty."); return

        print("Generating clustermap...")
        clustergrid = sns.clustermap(clustermap_matrix, cmap='vlag', center=0, figsize=(12, 14),
                                     annot=False, linewidths=0, xticklabels=False, yticklabels=False,
                                     dendrogram_ratio=(.1, .2), cbar_pos=(0.02, 0.8, .03, .15),
                                     cbar_kws={'label': f'Spearman {rho_col}'})
        clustergrid.fig.suptitle('Global Correlation Structure (Clustermap)', y=1.02)

        if save_path:
            filename = os.path.join(save_path, "global_correlation_clustermap.png")
            clustergrid.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved Clustermap to: {filename}")
        print("NOTE: Labels hidden due to size.")
        plt.show(); plt.close()

    except MemoryError: print("Error: MemoryError generating clustermap.")
    except Exception as e: print(f"Error generating clustermap: {e}"); plt.close()


def plot_top_correlation_heatmap(df_results, p_adj_col='p_adjusted', rho_col='rho',
                                 gene_col='gene', metab_col='metabolite',
                                 top_n_pairs=PLOT_HEATMAP_TOP_N_PAIRS,
                                 max_features=PLOT_HEATMAP_MAX_FEATURES,
                                 save_path=None):
    """Generates a heatmap for features involved in top correlations."""
    print("\n--- Generating Heatmap of Top Correlated Features ---")
    required_cols = [gene_col, metab_col, rho_col, p_adj_col]
    if df_results is None or df_results.empty or not all(c in df_results.columns for c in required_cols):
        print(f"Cannot plot: Input data invalid or missing columns: {required_cols}")
        return

    try:
        print(f"Selecting top {top_n_pairs} pairs by {p_adj_col}...")
        top_pairs = df_results.nsmallest(top_n_pairs, p_adj_col)
        top_gene_counts = top_pairs[gene_col].value_counts()
        top_metabolite_counts = top_pairs[metab_col].value_counts()
        top_genes = top_gene_counts.nlargest(max_features).index.tolist()
        top_metabolites = top_metabolite_counts.nlargest(max_features).index.tolist()
        print(f"Identified top {len(top_genes)} genes and {len(top_metabolites)} metabolites.")

        heatmap_data = df_results[df_results[gene_col].isin(top_genes) &
                                  df_results[metab_col].isin(top_metabolites)]
        if heatmap_data.empty: print("Skipping heatmap: No correlations between selected features."); return

        heatmap_pivot = heatmap_data.pivot_table(index=gene_col, columns=metab_col, values=rho_col)
        # Optional: Reindex for consistent ordering
        # heatmap_pivot = heatmap_pivot.reindex(index=top_genes, columns=top_metabolites)
        print(f"Plotting heatmap ({heatmap_pivot.shape[0]} x {heatmap_pivot.shape[1]})...")
        h_w = min(heatmap_pivot.shape[1] * 0.5 + 5, 30)
        h_h = min(heatmap_pivot.shape[0] * 0.4 + 3, 25)
        plt.figure(figsize=(h_w, h_h))
        sns.heatmap(heatmap_pivot, cmap='vlag', center=0, annot=False, linewidths=.5, cbar_kws={'label': f'Spearman {rho_col}'})
        plt.title(f'Heatmap of Top {max_features} Features (by frequency in top {top_n_pairs} pairs)')
        plt.xlabel('Metabolites', fontsize=10); plt.ylabel('Genes', fontsize=10)
        plt.xticks(rotation=90, fontsize=8); plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()

        if save_path:
            filename = os.path.join(save_path, f"top_{max_features}_features_correlations_heatmap.png")
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Saved Heatmap to: {filename}")
        plt.show(); plt.close()

    except Exception as e: print(f"Error generating Heatmap: {e}"); plt.close()


def plot_top_scatter(df_results, rna_data_filtered, metabolite_data_filtered,
                     top_n=PLOT_N_TOP_POS_NEG_SCATTER, rho_col='rho',
                     p_adj_col='p_adjusted', gene_col='gene', metab_col='metabolite',
                     save_path=None):
    """Generates scatter plots for top positive and top negative correlations."""
    print("\n--- Generating Top Correlation Scatter Plots ---")
    # Check prerequisites
    if df_results is None or df_results.empty: print("Cannot plot: Missing correlation results."); return
    if rna_data_filtered is None or rna_data_filtered.empty: print("Cannot plot: Missing filtered RNA data."); return
    if metabolite_data_filtered is None or metabolite_data_filtered.empty: print("Cannot plot: Missing filtered Metabolite data."); return
    if not rna_data_filtered.index.equals(metabolite_data_filtered.index): print("Cannot plot: Filtered data indices mismatch."); return
    required_cols = [gene_col, metab_col, rho_col, p_adj_col]
    if not all(c in df_results.columns for c in required_cols): print(f"Cannot plot: Results missing columns: {required_cols}"); return

    top_pos = df_results.nlargest(top_n, rho_col)
    top_neg = df_results.nsmallest(top_n, rho_col)

    for df_subset, label, color, filename_suffix in [
        (top_pos, "Positive", 'darkgreen', "positive"),
        (top_neg, "Negative", 'purple', "negative")
    ]:
        if df_subset.empty:
            print(f"Skipping {label.lower()} scatter plots: No correlations found.")
            continue

        print(f"Generating top {len(df_subset)} {label.lower()} correlation scatter plots...")
        try:
            num_plots = len(df_subset)
            ncols = min(num_plots, 5); nrows = (num_plots + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 4))
            axes = np.array(axes).flatten()

            for i, (idx, row) in enumerate(df_subset.iterrows()):
                gene, metab, rho, p_adj = row[gene_col], row[metab_col], row[rho_col], row[p_adj_col]
                ax = axes[i]
                if gene in rna_data_filtered.columns and metab in metabolite_data_filtered.columns:
                    sns.regplot(x=rna_data_filtered[gene], y=metabolite_data_filtered[metab], ax=ax,
                                scatter_kws={'alpha':0.7, 's': 30, 'color': color if label=='Positive' else None}, # Use specific color for pos, default (blue) for neg
                                line_kws={'color':'red', 'lw':1.5})
                    ax.set_title(f"{gene} vs {metab}\nRho={rho:.3f}, p.adj={p_adj:.2e}", fontsize=9)
                    ax.set_xlabel(f"{gene}", fontsize=8); ax.set_ylabel(f"{metab}", fontsize=8)
                    ax.tick_params(axis='both', which='major', labelsize=7); ax.grid(alpha=0.3)
                else:
                    ax.text(0.5, 0.5, f'Data missing', ha='center', va='center', color='red', fontsize=8)
                    ax.set_title(f"Data Missing", fontsize=9); ax.set_xticks([]); ax.set_yticks([])

            for j in range(num_plots, len(axes)): axes[j].set_visible(False)
            fig.suptitle(f"Top {len(df_subset)} {label} Correlations (by Rho)", y=1.03, fontsize=14)
            plt.tight_layout(rect=[0, 0, 1, 0.97])

            if save_path:
                plot_filename = os.path.join(save_path, f"top_{filename_suffix}_correlation_scatter_plots.png")
                plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
                print(f"Saved Top {label} scatter plots to: {plot_filename}")
            plt.show(); plt.close()

        except Exception as e: print(f"Error generating {label.lower()} scatter plots: {e}"); plt.close()


def evaluate_nmf(H_df, W_rna_df, W_metab_df, V_input_matrix, model, k, plots_save_path=None):
     """Calculates and prints NMF evaluation metrics."""
     print("\n--- Evaluating NMF Results ---")
     if H_df is None or W_rna_df is None or W_metab_df is None or model is None:
          print("Cannot evaluate: Missing NMF result components (H_df, W_rna_df, W_metab_df, model).")
          return None

     metrics = {}

     # Metric 1: Reconstruction Error
     print("\n  --- 1. Reconstruction Error ---")
     try:
        recon_error = model.reconstruction_err_
        metrics['Reconstruction Error'] = recon_error
        print(f"    - Frobenius Norm Error (from model): {recon_error:.4f}")
        if V_input_matrix is not None:
            norm_V = np.linalg.norm(V_input_matrix.values, 'fro')
            relative_error = recon_error / norm_V if norm_V > 0 else np.nan
            metrics['Relative Error'] = relative_error
            print(f"    - Relative Error (approx.): {relative_error:.4f}")
        else: print("    - Relative Error: Cannot calculate without input matrix V.")
     except Exception as e: print(f"Error calculating reconstruction error: {e}")

     # Metric 2: Component Contribution
     print("\n  --- 2. Component Contribution (W Matrix) ---")
     try:
        rna_loading_sq_sum = (W_rna_df**2).sum(axis=1); metab_loading_sq_sum = (W_metab_df**2).sum(axis=1)
        total_loading_sq_sum = rna_loading_sq_sum + metab_loading_sq_sum
        prop_rna = rna_loading_sq_sum / total_loading_sq_sum.replace(0, np.nan)
        prop_metab = metab_loading_sq_sum / total_loading_sq_sum.replace(0, np.nan)
        contribution_df = pd.DataFrame({'Prop_RNA': prop_rna,'Prop_Metab': prop_metab}, index=W_rna_df.index)
        print(contribution_df.to_string(float_format="%.4f"))
        metrics['Contribution DF'] = contribution_df
     except Exception as e: print(f"Error calculating contribution metrics: {e}")

     # Metric 3: Sparsity
     print("\n  --- 3. Sparsity ---")
     try:
        sparsity_threshold = 1e-8
        sparsity_H = np.mean(np.abs(H_df.values) < sparsity_threshold)
        sparsity_W_rna = np.mean(np.abs(W_rna_df.values) < sparsity_threshold)
        sparsity_W_metab = np.mean(np.abs(W_metab_df.values) < sparsity_threshold)
        metrics['Sparsity H'] = sparsity_H; metrics['Sparsity W RNA'] = sparsity_W_rna; metrics['Sparsity W Metab'] = sparsity_W_metab
        print(f"    - Sparsity H: {sparsity_H:.4f} | W_rna: {sparsity_W_rna:.4f} | W_metab: {sparsity_W_metab:.4f}")
     except Exception as e: print(f"Error calculating sparsity metrics: {e}")

     # Metric 4: H Correlation
     print("\n  --- 4. H Matrix Correlation ---")
     try:
        if H_df.shape[1] > 1:
            corr_H = H_df.corr()
            np.fill_diagonal(corr_H.values, np.nan)
            max_h_corr = corr_H.abs().max().max()
            metrics['Max H Correlation'] = max_h_corr
            print(f"    - Max absolute correlation between components: {max_h_corr:.4f}")
            if k <= 15:
                 plt.figure(figsize=(6, 5))
                 sns.heatmap(corr_H, cmap='coolwarm', center=0, annot=True, fmt=".2f", linewidths=.5)
                 plt.title(f'Correlation of Sample Loadings (H, k={k})')
                 if plots_save_path:
                      plot_filename_h_corr = os.path.join(plots_save_path, f"nmf_H_correlation_k{k}.png")
                      try: plt.savefig(plot_filename_h_corr, dpi=150, bbox_inches='tight'); print(f"    Saved H correlation plot to: {plot_filename_h_corr}")
                      except Exception as e_p: print(f"    Error saving plot: {e_p}")
                 plt.show(); plt.close()
            else: print("    (Skipping H correlation heatmap, k > 15)")
        else: print("    - Skipping correlation (k=1).")
     except Exception as e: print(f"Error calculating component correlation: {e}")

     return metrics