# mercat_analyzer/mercat_analyzer/analysis.py

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.decomposition import NMF
from sklearn.preprocessing import MaxAbsScaler
import time
from .config import NMF_DEFAULT_K, NMF_DEFAULT_MAX_ITER, NMF_DEFAULT_RANDOM_STATE

def run_spearman_correlation(rna_data_filtered, metabolite_data_filtered):
    """Calculates pairwise Spearman correlations between two aligned dataframes."""
    print("\n--- Running Spearman Correlation Analysis ---")
    if rna_data_filtered is None or metabolite_data_filtered is None or \
       rna_data_filtered.empty or metabolite_data_filtered.empty:
        print("Error: Input filtered RNA or Metabolite data is missing or empty.")
        return None
    if not rna_data_filtered.index.equals(metabolite_data_filtered.index):
         print("Error: Indices of input dataframes do not match!")
         return None

    print(f"Input RNA shape: {rna_data_filtered.shape}")
    print(f"Input Metabolite shape: {metabolite_data_filtered.shape}")

    try:
        # Ensure data is numeric for correlation
        rna_numeric = rna_data_filtered.select_dtypes(include=np.number)
        metab_numeric = metabolite_data_filtered.select_dtypes(include=np.number)
        if rna_numeric.empty or metab_numeric.empty:
             raise ValueError("No numeric columns found in one or both inputs for correlation.")

        metabolite_array = metab_numeric.values
        rna_array = rna_numeric.values
        print("Converted data to numpy arrays.")

        metabolite_names = metab_numeric.columns.tolist()
        gene_names = rna_numeric.columns.tolist()
        n_genes = rna_array.shape[1]
        n_mets = metabolite_array.shape[1]

        print(f"Starting correlations: {n_genes} genes x {n_mets} metabolites = {n_genes * n_mets} pairs.")
        start_time = time.time()

        results = [] # Store tuples (gene, metabolite, rho, p)
        for j in range(n_mets):
            metabolite_vector = metabolite_array[:, j]
            for i in range(n_genes):
                try:
                    rho, p = spearmanr(rna_array[:, i], metabolite_vector)
                except ValueError: # Handle cases like zero variance
                    rho, p = np.nan, np.nan
                results.append((gene_names[i], metabolite_names[j], rho, p))

        duration = time.time() - start_time
        print(f"Correlation calculations finished in {duration:.2f} seconds.")

        # Create DataFrame
        df_res = pd.DataFrame(results, columns=['gene', 'metabolite', 'rho', 'p'])
        df_res.dropna(subset=['p', 'rho'], inplace=True) # Remove failed calculations
        print(f"Generated {len(df_res)} valid correlation pairs.")

        return df_res

    except Exception as e:
        print(f"An error occurred during Spearman correlation: {e}")
        return None


def run_nmf_concatenated(rna_data_filtered, metabolite_data_filtered,
                         n_components=NMF_DEFAULT_K,
                         max_iter=NMF_DEFAULT_MAX_ITER,
                         random_state=NMF_DEFAULT_RANDOM_STATE):
    """Runs concatenated NMF (Strategy 2) on aligned, filtered, non-negative data."""
    print("\n--- Running Concatenated NMF Analysis ---")
    nmf_results = {
        "H_df": None, "W_rna_df": None, "W_metab_df": None,
        "model": None, "V_combined_imputed": None, "successful": False
    }

    # --- 1. Prerequisites & Non-Negativity ---
    prerequisites_met = True
    if rna_data_filtered is None or metabolite_data_filtered is None or \
       rna_data_filtered.empty or metabolite_data_filtered.empty:
        print("ERROR: Input filtered data missing or empty."); prerequisites_met = False
    elif not rna_data_filtered.index.equals(metabolite_data_filtered.index):
        print("ERROR: Input indices do not match!"); prerequisites_met = False

    rna_numeric = None; metab_numeric = None
    if prerequisites_met:
        try:
            rna_numeric = rna_data_filtered.select_dtypes(include=np.number)
            metab_numeric = metabolite_data_filtered.select_dtypes(include=np.number)
            if rna_numeric.empty or metab_numeric.empty: raise ValueError("No numeric columns.")
            if rna_numeric.min().min() < -1e-9 or metab_numeric.min().min() < -1e-9: raise ValueError("Negative values detected.")
            print("Prerequisites and non-negativity checks passed.")
            original_rna_features = rna_numeric.columns.tolist()
            original_metab_features = metab_numeric.columns.tolist()
            n_rna_features = len(original_rna_features)
            n_metab_features = len(original_metab_features)
            sample_index = rna_numeric.index
        except Exception as e: print(f"Error in prerequisite check: {e}"); prerequisites_met = False

    # --- 2. Scaling ---
    rna_scaled_df = None; metab_scaled_df = None
    if prerequisites_met:
        print("Scaling data (MaxAbsScaler)...")
        try:
            scaler_rna = MaxAbsScaler(); rna_scaled_array = scaler_rna.fit_transform(rna_numeric.values)
            rna_scaled_df = pd.DataFrame(rna_scaled_array, index=sample_index, columns=original_rna_features)
            scaler_metab = MaxAbsScaler(); metab_scaled_array = scaler_metab.fit_transform(metab_numeric.values)
            metab_scaled_df = pd.DataFrame(metab_scaled_array, index=sample_index, columns=original_metab_features)
            if rna_scaled_df.min().min() < -1e-9 or metab_scaled_df.min().min() < -1e-9: print("Warning: Negatives after scaling?")
        except Exception as e: print(f"Error during scaling: {e}"); prerequisites_met = False

    # --- 3. NaN Handling & Concatenation ---
    V_combined_imputed = None
    if prerequisites_met:
        print("Handling NaNs and concatenating...")
        try:
            V_combined_scaled = pd.concat([rna_scaled_df, metab_scaled_df], axis=1)
            nan_count = V_combined_scaled.isna().sum().sum()
            if nan_count > 0:
                print(f"Imputing {nan_count} NaNs with 0...")
                V_combined_imputed = V_combined_scaled.fillna(0)
                if V_combined_imputed.isna().sum().sum() != 0: raise ValueError("NaN imputation failed!")
            else: V_combined_imputed = V_combined_scaled; print("No NaNs found.")
            if V_combined_imputed.min().min() < -1e-9: raise ValueError("Negative values after imputation!")
            print(f"Final matrix shape for NMF: {V_combined_imputed.shape}")
        except Exception as e: print(f"Error during NaN/Concat: {e}"); prerequisites_met = False

    # --- 4. Run NMF ---
    if prerequisites_met:
        print(f"Running NMF (k={n_components}, max_iter={max_iter})...")
        nmf_start_time = time.time()
        try:
            model = NMF(n_components=n_components, init='nndsvda', max_iter=max_iter,
                        random_state=random_state, solver='cd', beta_loss='frobenius', tol=1e-4)
            H_combined = model.fit_transform(V_combined_imputed.values)
            W_combined = model.components_
            duration = time.time() - nmf_start_time
            print(f"NMF fitting completed in {duration:.2f}s. Recon Error: {model.reconstruction_err_:.4f}")

            # Process output
            component_names = [f"Comp_{i+1}" for i in range(n_components)]
            H_df = pd.DataFrame(H_combined, index=sample_index, columns=component_names)
            W_combined_df = pd.DataFrame(W_combined, index=component_names, columns=V_combined_imputed.columns)
            W_rna_df = W_combined_df[original_rna_features].copy()
            W_metab_df = W_combined_df[original_metab_features].copy()

            # Store results
            nmf_results["H_df"] = H_df
            nmf_results["W_rna_df"] = W_rna_df
            nmf_results["W_metab_df"] = W_metab_df
            nmf_results["model"] = model
            nmf_results["V_combined_imputed"] = V_combined_imputed # Save the input matrix too
            nmf_results["successful"] = True
            print("NMF results processed.")

        except Exception as e_nmf: print(f"Error during NMF execution: {e_nmf}")

    if not nmf_results["successful"]: print("NMF run failed or prerequisites not met.")
    return nmf_results