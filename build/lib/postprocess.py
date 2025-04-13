# mercat_analyzer/mercat_analyzer/postprocess.py

import pandas as pd
import numpy as np
try:
    from statsmodels.sandbox.stats.multicomp import multipletests
    statsmodels_available = True
except ImportError:
    statsmodels_available = False

def adjust_pvalues_bh(df_results, p_col='p', q_col='p_adjusted'):
    """Applies Benjamini-Hochberg FDR correction."""
    print("\n--- Applying Benjamini-Hochberg FDR Correction ---")
    if df_results is None or df_results.empty:
        print("Input results DataFrame is missing or empty. Skipping.")
        return None
    if p_col not in df_results.columns:
        print(f"Error: Raw p-value column '{p_col}' not found.")
        return df_results # Return original df

    if not statsmodels_available:
        print("Warning: 'statsmodels' not found. Cannot perform correction.")
        return df_results

    try:
        df_adj = df_results.copy()
        print(f"Processing {len(df_adj)} p-values...")

        # Handle NaNs and ensure range [0, 1]
        p_values_clean = df_adj[p_col].fillna(1.0).clip(0.0, 1.0)
        nan_count = df_adj[p_col].isna().sum()
        if nan_count > 0: print(f"  Filled {nan_count} NaNs with 1.0.")

        # Perform correction
        reject, pvals_corrected, _, _ = multipletests(p_values_clean, method='fdr_bh')

        df_adj[q_col] = pvals_corrected
        df_adj = df_adj.sort_values(by=q_col, ascending=True)
        print(f"Added '{q_col}' column and sorted results.")
        return df_adj

    except Exception as e:
        print(f"Error during p-value adjustment: {e}")
        # Return original dataframe if adjustment fails
        return df_results