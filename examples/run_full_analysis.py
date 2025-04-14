# Example Script using the MeERCAT package (Local File Version - Final Alignment Fix)

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# --- Add Package to Path (If necessary) ---
try:
    package_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
    if package_path not in sys.path: sys.path.insert(0, package_path)
except NameError: print("Assuming 'meercat' is installed or accessible.")

# Import functions from YOUR package
try:
    import meercat; from meercat import config, utils, load_data, preprocess, analysis, postprocess, visualize
except ImportError as e: print(f"ERROR importing 'meercat': {e}"); sys.exit(1)

# ==================================================
# --- User Configuration ---
# ==================================================
BASE_PROJECT_PATH = '.'
INPUT_DATA_SUBDIR = 'input_data'
OUTPUT_RESULTS_SUBDIR = 'output_results'
# ==================================================


# --- 1. Setup Paths ---
print("--- 1. Setting up Paths ---")
input_data_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, INPUT_DATA_SUBDIR))
output_results_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, OUTPUT_RESULTS_SUBDIR))
paths = utils.setup_paths(output_results_path, config.DEFAULT_PROJECT_FOLDER, config.DEFAULT_SPEARMAN_SUBFOLDER, config.DEFAULT_NMF_SUBFOLDER, config.DEFAULT_PLOTS_SUBFOLDER)
paths['input'] = input_data_path; paths['input_metabolites'] = os.path.join(input_data_path, 'metabolites')
paths['input_rna'] = os.path.join(input_data_path, 'rna'); paths['input_rna_metadata'] = os.path.join(input_data_path, config.RNA_METADATA_FILENAME)
print(f"\nInput Data expected in: {paths['input']}") # etc...
print(f"Output Results will be saved under: {paths['base']}")
save_files = True


# --- 2. Load Data ---
print("\n--- 2. Loading Data ---")
loaded_metabolite_data = load_data.load_metabolite_files(paths['input_metabolites'])
loaded_rna_data = load_data.load_rna_files(paths['input_rna'], rows_are_genes=True)
external_rna_metadata_df = load_data.load_external_rna_metadata(paths['input_rna_metadata'])
# --- Exit if critical data failed ---
critical_data_missing = False
if not loaded_metabolite_data: print("CRITICAL ERROR: No metabolite data."); critical_data_missing = True
if not loaded_rna_data: print("CRITICAL ERROR: No RNA counts."); critical_data_missing = True
if external_rna_metadata_df is None: print("CRITICAL ERROR: External RNA metadata failed load."); critical_data_missing = True
if critical_data_missing: print("Exiting."); sys.exit(1)


# --- 3. Preprocess Metabolites ---
print("\n--- 3. Preprocessing Metabolites ---")
metabolite_features = None # Features indexed by original_rna_sample_id
metabolite_metadata = None # Metadata indexed by original_rna_sample_id
metabolite_combined_raw = preprocess.combine_dataframes(loaded_metabolite_data, axis=0)
if metabolite_combined_raw is not None:
    # clean_metabolite_data returns features and metadata indexed by original_rna_sample_id (if unique)
    metabolite_features, metabolite_metadata = preprocess.clean_metabolite_data(metabolite_combined_raw)
    # *** Check index name AFTER calling clean_metabolite_data ***
    if metabolite_features is None or metabolite_features.index.name != 'original_rna_sample_id':
         print("ERROR: Metabolite preprocessing failed or index not set to 'original_rna_sample_id'.")
         metabolite_features = None # Invalidate if index is wrong
    else:
         print("Metabolite preprocessing complete (indexed by original_rna_sample_id).")
         if save_files: # Save if successful
              metab_feat_path = os.path.join(paths['base'], "metabolite_features_originalID.csv") # Adjusted name
              meta_extract_path = os.path.join(paths['base'], "metabolite_metadata_originalID.csv") # Adjusted name
              try: metabolite_features.to_csv(metab_feat_path, index=True); print(f"Saved metabolite FEATURES to {metab_feat_path}")
              except Exception as e: print(f"Warning: Could not save metabolite features: {e}")
              if metabolite_metadata is not None: # Metadata might exist even if index failed on features
                  try: metabolite_metadata.to_csv(meta_extract_path, index=True); print(f"Saved metabolite METADATA to {meta_extract_path}")
                  except Exception as e: print(f"Warning: Could not save metabolite metadata: {e}")
else: print("Skipping metabolite preprocessing - no combined raw data.")


# --- 4. Preprocess RNA ---
print("\n--- 4. Preprocessing RNA ---")
rna_normalized = None      # Final normalized features indexed by original_rna_sample_id
rna_combined_features = None # Combined raw features indexed by original_rna_sample_id
processed_rna_features_dict = {}

if loaded_rna_data:
    for name, df_raw in loaded_rna_data.items():
        # process_rna_dataframe returns features indexed by original sample ID
        features_df = preprocess.process_rna_dataframe(df_raw, name, rows_are_genes=True)
        if features_df is not None: processed_rna_features_dict[name] = features_df

# Combine FEATURES ONLY (indexed by original sample ID)
rna_combined_features = preprocess.combine_dataframes(processed_rna_features_dict, axis=0, join='outer')

# Check overlap with External Metadata (sanity check)
if rna_combined_features is not None and external_rna_metadata_df is not None:
    print("\nChecking overlap between RNA features and external metadata...")
    common_rna_ids = rna_combined_features.index.intersection(external_rna_metadata_df.index)
    print(f"Found {len(common_rna_ids)} matching sample IDs based on original_rna_sample_id.")
    if len(common_rna_ids) == 0:
        print("ERROR: No overlap between RNA features and metadata. Check input files.")
        rna_combined_features = None # Invalidate if no overlap
    elif len(common_rna_ids) < len(rna_combined_features.index):
        print("WARNING: Subsetting RNA features to match metadata.")
        rna_combined_features = rna_combined_features.loc[common_rna_ids] # Subset features

# Normalize RNA Data (use features indexed by original_rna_sample_id)
if rna_combined_features is not None:
    # normalize_rna_data expects input indexed by original_rna_sample_id
    if rna_combined_features.index.name != 'original_rna_sample_id':
         print("ERROR: Combined RNA features index name is incorrect before normalization!")
         rna_normalized = None
    else:
        rna_normalized = preprocess.normalize_rna_data(rna_combined_features) # Normalize features
        if rna_normalized is not None:
            print(f"RNA Normalization complete. Shape: {rna_normalized.shape}")
            if save_files:
                 rna_norm_path = os.path.join(paths['base'], config.RNA_NORMALIZED_FILENAME) # Save normalized features
                 try: rna_normalized.to_csv(rna_norm_path, index=True); print(f"Saved normalized RNA FEATURES to {rna_norm_path}")
                 except Exception as e: print(f"Warning: Could not save normalized RNA data: {e}")
        else: print("RNA normalization failed.")
else: print("Skipping RNA normalization - combined RNA features invalid or missing.")


print("\n--- 5. Aligning Samples ---")
metabolite_matched = None; rna_matched = None

# *** Perform FINAL cleaning and type check right before alignment ***
print("\n--- DEBUG: Final Index Cleaning & Type Check Before Alignment ---")
metab_ok_final = False
rna_ok_final = False
temp_metab_features = None
temp_rna_normalized = None

if 'metabolite_features' in locals() and metabolite_features is not None and not metabolite_features.empty:
    if metabolite_features.index.name == 'original_rna_sample_id':
        try:
            # Create temporary copies for final cleaning
            temp_metab_features = metabolite_features.copy()
            # Explicitly clean the index AGAIN
            temp_metab_features.index = preprocess.clean_index(temp_metab_features.index)
            temp_metab_features.index.name = 'original_rna_sample_id' # Re-assign name
            print("Metabolite Index FINAL Cleaned (Top 10):", temp_metab_features.index[:10].tolist())
            print(f"Metabolite Index FINAL dtype: {temp_metab_features.index.dtype}")
            metab_ok_final = True
        except Exception as e_clean_m:
            print(f"ERROR during final metabolite index cleaning: {e_clean_m}")
    else: print("ERROR: Metabolite features FINAL missing correct index name.")
else: print("Metabolite Features FINAL missing or empty.")

print("-" * 20)

if 'rna_normalized' in locals() and rna_normalized is not None and not rna_normalized.empty:
     if rna_normalized.index.name == 'original_rna_sample_id':
        try:
            # Create temporary copies for final cleaning
            temp_rna_normalized = rna_normalized.copy()
            # Explicitly clean the index AGAIN
            temp_rna_normalized.index = preprocess.clean_index(temp_rna_normalized.index)
            temp_rna_normalized.index.name = 'original_rna_sample_id' # Re-assign name
            print("RNA Normalized Index FINAL Cleaned (Top 10):", temp_rna_normalized.index[:10].tolist())
            print(f"RNA Normalized Index FINAL dtype: {temp_rna_normalized.index.dtype}")
            rna_ok_final = True
        except Exception as e_clean_r:
             print(f"ERROR during final RNA index cleaning: {e_clean_r}")
     else: print("ERROR: RNA Normalized FINAL missing correct index name.")
else: print("RNA Normalized FINAL missing or empty.")

print("--- END FINAL DEBUG ---")


# *** Call align_samples using the TEMP cleaned dataframes ***
if metab_ok_final and rna_ok_final:
    metabolite_matched, rna_matched = preprocess.align_samples(
        temp_metab_features, # Use the re-cleaned version
        temp_rna_normalized, # Use the re-cleaned version
        df1_name="Metabolites", df2_name="RNA"
    )
    # ... (rest of alignment saving/checks using metabolite_matched, rna_matched) ...
    if metabolite_matched is not None and rna_matched is not None:
         print(f"Alignment successful. Matched {len(metabolite_matched)} samples.")
         # ... (saving) ...
    else: print("ERROR: Sample alignment function returned None."); sys.exit(1)
else:
    print("Skipping Alignment: Final check/cleaning failed for metabolite or RNA features.")
    sys.exit(1)



    
# --- 6. Variance Filtering ---
print("\n--- 6. Variance Filtering ---")
# Inputs are metabolite_matched, rna_matched (indexed by original_rna_sample_id)
rna_data_filtered = preprocess.apply_variance_filter(rna_matched, config.VAR_FILTER_TOP_N_GENES, "RNA")
metabolite_data_filtered = preprocess.apply_variance_filter(metabolite_matched, config.VAR_FILTER_TOP_N_METABOLITES, "Metabolite")
if rna_data_filtered is None or metabolite_data_filtered is None: print("ERROR: Variance filtering failed."); sys.exit(1)
# NOTE: Output dataframes are still indexed by original_rna_sample_id


# --- 7. Spearman Correlation ---
print("\n--- 7. Spearman Correlation ---")
# Inputs are rna_data_filtered, metabolite_data_filtered
df_corr_raw = analysis.run_spearman_correlation(rna_data_filtered, metabolite_data_filtered)
if df_corr_raw is None: print("ERROR: Spearman correlation failed."); sys.exit(1)
if save_files: # Save raw correlations
     corr_raw_path = os.path.join(paths['spearman'], config.CORRELATION_RAW_FILENAME)
     try: df_corr_raw.to_csv(corr_raw_path, index=False); print(f"Saved raw Spearman results to {corr_raw_path}")
     except Exception as e: print(f"Warning: Could not save raw correlation results: {e}")


# --- 8. P-value Adjustment ---
print("\n--- 8. P-value Adjustment ---")
df_corr_adj = postprocess.adjust_pvalues_bh(df_corr_raw)
if df_corr_adj is None: print("ERROR: P-value adjustment failed."); df_corr_adj = df_corr_raw
if save_files: # Save adjusted correlations
     corr_adj_path = os.path.join(paths['spearman'], config.CORRELATION_ADJ_FILENAME)
     try: df_corr_adj.to_csv(corr_adj_path, index=False); print(f"Saved adjusted Spearman results to {corr_adj_path}")
     except Exception as e: print(f"Warning: Could not save adjusted correlation results: {e}")


# --- 9. Correlation Visualization ---
print("\n--- 9. Correlation Visualization ---")
if df_corr_adj is not None:
     print(f"Generating correlation plots (saving to {paths['spearman_plots']})...")
     visualize.plot_rho_distribution(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_volcano(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_correlation_clustermap(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_top_correlation_heatmap(df_corr_adj, save_path=paths['spearman_plots'])
     # Pass the variance filtered data (still indexed by original_rna_sample_id)
     visualize.plot_top_scatter(df_corr_adj, rna_data_filtered, metabolite_data_filtered, save_path=paths['spearman_plots'])
else: print("Skipping correlation visualization.")


# --- 10. NMF Analysis & Evaluation ---
print("\n--- 10. NMF Analysis & Evaluation ---")
k_values_to_test = [3, 4, 5, 6, 7]
all_nmf_metrics = {}
for k_to_run in k_values_to_test:
    print("\n" + "*"*30); print(f" Running NMF for k = {k_to_run}"); print("*"*30)
    if rna_data_filtered is None or metabolite_data_filtered is None: print(f"Skipping NMF k={k_to_run}: Filtered data missing."); continue
    # NMF input data is indexed by original_rna_sample_id
    nmf_results = analysis.run_nmf_concatenated(
        rna_data_filtered, metabolite_data_filtered, n_components=k_to_run,
        max_iter=config.NMF_DEFAULT_MAX_ITER, random_state=config.NMF_DEFAULT_RANDOM_STATE
    )
    if nmf_results["successful"]:
        # NMF H matrix index will be original_rna_sample_id
        if save_files: # Save results
             try:
                 h_filename=config.NMF_H_FILENAME_TEMPLATE.format(k=k_to_run); w_rna_filename=config.NMF_W_RNA_FILENAME_TEMPLATE.format(k=k_to_run); w_metab_filename=config.NMF_W_METAB_FILENAME_TEMPLATE.format(k=k_to_run)
                 if nmf_results["H_df"] is not None: nmf_results["H_df"].to_csv(os.path.join(paths['nmf'], h_filename), index=True)
                 if nmf_results["W_rna_df"] is not None: nmf_results["W_rna_df"].to_csv(os.path.join(paths['nmf'], w_rna_filename), index=True)
                 if nmf_results["W_metab_df"] is not None: nmf_results["W_metab_df"].to_csv(os.path.join(paths['nmf'], w_metab_filename), index=True)
                 print(f"Saved NMF results for k={k_to_run} to {paths['nmf']}")
             except Exception as e: print(f"Warning: Could not save NMF results k={k_to_run}: {e}")
        print(f"\n--- Evaluating NMF k={k_to_run} ---") # Evaluate results
        # Pass variance filtered data for V recalc (indexed by original_rna_sample_id)
        evaluation_metrics = visualize.evaluate_nmf(
            nmf_results["H_df"], nmf_results["W_rna_df"], nmf_results["W_metab_df"],
            rna_data_filtered=rna_data_filtered, metabolite_data_filtered=metabolite_data_filtered,
            model=nmf_results["model"], k=k_to_run, plots_save_path=paths['nmf_plots']
        )
        if evaluation_metrics: all_nmf_metrics[k_to_run] = evaluation_metrics; print(f"Stored evaluation metrics for k={k_to_run}")
    else: print(f"NMF analysis failed for k={k_to_run}.")

# --- Post-NMF Evaluation Plot ---
# (Keep plotting logic as before, including initializations)
k_vals = []; recon_errors = []; relative_errors = []
# ... (rest of post-NMF plotting) ...

print("\n--- FULL ANALYSIS PIPELINE COMPLETE ---")