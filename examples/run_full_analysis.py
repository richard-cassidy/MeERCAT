# Example Script using the MeERCAT package (Local File Version - Updated for Preprocess Refactor)

import pandas as pd
import numpy as np
import os
import sys # To potentially add package path if not installed
import matplotlib.pyplot as plt # Often needed to manage plots

# --- Add Package to Path (If necessary) ---
# If you haven't formally installed the package using 'pip install .',
# you might need to tell Python where to find it.
# Adjust the path '../' if your script is not directly inside the 'examples' folder.
try:
    # Assumes the script is in a subdirectory (like 'examples') relative to the package root
    package_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
    if package_path not in sys.path:
        print(f"Adding package path: {package_path}")
        sys.path.insert(0, package_path)
except NameError:
     # __file__ is not defined in interactive environments like basic Python shell or some IDEs
     # Assume package is installed or in current working directory structure
     print("Running in environment where __file__ is not defined. Assuming 'meercat' is installed or accessible.")


# Import functions from YOUR package
try:
    import meercat
    from meercat import config
    from meercat import utils
    from meercat import load_data
    from meercat import preprocess # This now contains the refactored functions
    from meercat import analysis
    from meercat import postprocess
    from meercat import visualize
except ImportError as e:
     print(f"ERROR: Could not import the 'meercat' package.")
     print(f"Ensure the package is installed (`pip install .` in the root directory)")
     print(f"or the path is correctly added to sys.path if running from source.")
     print(f"Import error details: {e}")
     sys.exit(1) # Exit if package not found

# ==================================================
# --- User Configuration ---
# ==================================================
# --- Define Base Path for Input Data and Output Results ---
# --- *** USER MUST SET THIS PATH *** ---
BASE_PROJECT_PATH = '.' # Use '.' for the current directory

# --- Define Subdirectory Names ---
INPUT_DATA_SUBDIR = 'input_data'   # Expects subfolders 'metabolites', 'rna', and 'rna_metadata.csv' inside this
OUTPUT_RESULTS_SUBDIR = 'output_results' # Where all analysis results will be saved
# ==================================================


# --- 1. Setup Paths ---
print("--- 1. Setting up Paths ---")
# Construct full paths based on user config
input_data_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, INPUT_DATA_SUBDIR))
output_results_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, OUTPUT_RESULTS_SUBDIR))

# Use the setup_paths utility from the package to create output structure
paths = utils.setup_paths(
    output_results_path, # Base path for *outputs*
    '',                  # Project folder name within output (use empty)
    config.DEFAULT_SPEARMAN_SUBFOLDER,
    config.DEFAULT_NMF_SUBFOLDER,
    config.DEFAULT_PLOTS_SUBFOLDER # This is relative within spearman/nmf folders
)
# Add specific input paths for clarity
paths['input'] = input_data_path
paths['input_metabolites'] = os.path.join(input_data_path, 'metabolites')
paths['input_rna'] = os.path.join(input_data_path, 'rna')
# *** Define path for the NEW required RNA metadata file ***
paths['input_rna_metadata'] = os.path.join(input_data_path, 'rna_metadata.csv') # Needs to exist!
# Optional main metadata file (not currently used in this script)
# paths['input_metadata'] = os.path.join(input_data_path, 'metadata.csv')

print(f"\nInput Data expected in: {paths['input']}")
print(f"  (Metabolites in: {paths['input_metabolites']})")
print(f"  (RNA counts in: {paths['input_rna']})")
print(f"  (RNA metadata in: {paths['input_rna_metadata']})")
print(f"Output Results will be saved under: {output_results_path}")
print(f"  Spearman Results: {paths['spearman']}")
print(f"  NMF Results: {paths['nmf']}")
print(f"  Spearman Plots: {paths['spearman_plots']}")
print(f"  NMF Eval Plots: {paths['nmf_plots']}")

# Flag for saving intermediate/final files
save_files = True # Set to False to only run analysis in memory


# --- 2. Load Data ---
print("\n--- 2. Loading Data ---")
# --- Load Metabolites ---
if not os.path.isdir(paths['input_metabolites']):
     print(f"ERROR: Metabolite input directory not found: {paths['input_metabolites']}")
     loaded_metabolite_data = {}
else:
     metabolite_file_dict = {f: os.path.join(paths['input_metabolites'], f)
                              for f in os.listdir(paths['input_metabolites']) if f.lower().endswith('.csv')}
     loaded_metabolite_data = load_data.load_metabolite_files(metabolite_file_dict)
     if not loaded_metabolite_data: print("ERROR: No metabolite files loaded.")

# --- Load RNA Counts ---
if not os.path.isdir(paths['input_rna']):
     print(f"ERROR: RNA input directory not found: {paths['input_rna']}")
     loaded_rna_data = {}
else:
     rna_file_dict = {f: os.path.join(paths['input_rna'], f)
                       for f in os.listdir(paths['input_rna']) if f.lower().endswith('.csv')}
     # *** IMPORTANT: Set rows_are_genes correctly for your data format ***
     loaded_rna_data = load_data.load_rna_files(rna_file_dict, rows_are_genes=True)
     if not loaded_rna_data: print("ERROR: No RNA files loaded.")

# --- Load External RNA Metadata ---
external_rna_metadata_df = load_data.load_external_rna_metadata(paths['input_rna_metadata'])
if external_rna_metadata_df is None:
    print("ERROR: Cannot proceed without external RNA metadata file. Please create 'rna_metadata.csv'.")
    sys.exit(1) # Exit if essential metadata is missing


# --- 3. Preprocess Metabolites ---
print("\n--- 3. Preprocessing Metabolites ---")
metabolite_features = None # Will hold features-only data with composite index
metabolite_metadata = None # Will hold metadata with composite index

metabolite_combined_raw = preprocess.combine_dataframes(loaded_metabolite_data, axis=0)
if metabolite_combined_raw is not None:
    # *** UPDATED: Capture both returned dataframes ***
    metabolite_features, metabolite_metadata = preprocess.clean_metabolite_data(metabolite_combined_raw)

    # Save cleaned data (optional intermediate step)
    if metabolite_features is not None and save_files:
        metab_feat_path = os.path.join(paths['base'], config.METABOLITE_CLEANED_FILENAME)
        try:
            metabolite_features.to_csv(metab_feat_path, index=True)
            print(f"Saved cleaned metabolite FEATURES to {metab_feat_path}")
        except Exception as e: print(f"Warning: Could not save cleaned metabolite features: {e}")
    if metabolite_metadata is not None and save_files:
        meta_extract_path = os.path.join(paths['base'], config.METADATA_EXTRACTED_FILENAME)
        try:
            metabolite_metadata.to_csv(meta_extract_path, index=True)
            print(f"Saved extracted metabolite METADATA to {meta_extract_path}")
        except Exception as e: print(f"Warning: Could not save extracted metabolite metadata: {e}")
else:
    print("Skipping metabolite preprocessing - no combined raw data.")


# --- 4. Preprocess RNA ---
print("\n--- 4. Preprocessing RNA ---")
rna_normalized = None      # Will hold final normalized features with composite index
rna_metadata_final = None  # Will hold final RNA metadata with composite index
rna_features_indexed = None # Intermediate step: Features with composite index

processed_rna_features_dict = {} # Store features only (indexed by original ID)

if loaded_rna_data: # Check if dict has content
    for name, df_raw in loaded_rna_data.items():
        # process_rna_dataframe now ONLY returns features indexed by original sample ID
        features_df = preprocess.process_rna_dataframe(df_raw, name, rows_are_genes=True) # Set rows_are_genes correctly
        if features_df is not None:
            processed_rna_features_dict[name] = features_df

# Combine FEATURES ONLY
rna_combined_features = preprocess.combine_dataframes(processed_rna_features_dict, axis=0, join='outer')

# Merge Features with External Metadata
rna_merged_data = None # Will hold features + metadata, indexed by original ID
if rna_combined_features is not None and external_rna_metadata_df is not None:
    print("\nMerging combined RNA features with external metadata...")
    try:
        # Ensure indices are aligned before merge
        common_rna_ids = rna_combined_features.index.intersection(external_rna_metadata_df.index)
        if len(common_rna_ids) == 0: raise ValueError("No common IDs between combined RNA features and external metadata!")
        print(f"Found {len(common_rna_ids)} matching RNA sample IDs for metadata merge.")

        rna_combined_features_common = rna_combined_features.loc[common_rna_ids]
        external_rna_metadata_common = external_rna_metadata_df.reindex(common_rna_ids) # Align metadata index

        rna_merged_data = pd.merge(rna_combined_features_common, external_rna_metadata_common,
                                   left_index=True, right_index=True, how='left')
        print(f"Shape after merging RNA features and metadata: {rna_merged_data.shape}")

    except Exception as e_merge:
        print(f"Error merging RNA features and metadata: {e_merge}")

# Create Composite Index on Merged RNA Data
if rna_merged_data is not None:
    print("\nCreating composite index for merged RNA data...")
    # Use the _create_composite_id helper (it's internal to preprocess but called via clean_metabolite_data)
    # We need similar logic here or call the helper explicitly if made public/refactored
    # Replicating logic for now:
    composite_id_series_rna = preprocess._create_composite_id(rna_merged_data, config.METADATA_COMPOSITE_ID_COLS)

    if composite_id_series_rna is not None:
         if composite_id_series_rna.duplicated().any():
             print("CRITICAL WARNING: Duplicate composite IDs generated for RNA AFTER merge! Check metadata file content.")
         else:
             print("RNA composite ID unique. Setting index and separating features/metadata...")
             # Set index on the merged data first
             rna_merged_data.index = preprocess.clean_index(composite_id_series_rna)
             rna_merged_data.index.name = 'composite_sample_id'
             # Separate features using known feature names (important!)
             feature_cols_rna = rna_combined_features.columns # Use columns from before merge
             rna_features_indexed = rna_merged_data[feature_cols_rna].copy()
             # Separate metadata
             metadata_cols_rna = external_rna_metadata_df.columns.tolist() # Get original metadata cols
             # Add experiment_id if it wasn't in external file but was generated by file loading
             if 'experiment_id' in rna_merged_data.columns and 'experiment_id' not in metadata_cols_rna:
                  metadata_cols_rna.append('experiment_id')
             # Get metadata columns actually present in merged data after indexing
             final_rna_meta_cols = [col for col in metadata_cols_rna if col in rna_merged_data.columns]
             rna_metadata_final = rna_merged_data[final_rna_meta_cols].copy()
             print(f"Set composite index. RNA Features shape: {rna_features_indexed.shape}, RNA Metadata shape: {rna_metadata_final.shape}")
    else:
         print("Failed to create composite ID for merged RNA data.")

# Normalize RNA Data (use features with composite index)
if rna_features_indexed is not None:
    rna_normalized = preprocess.normalize_rna_data(rna_features_indexed) # Normalize indexed features
    # Save normalized data
    if rna_normalized is not None and save_files:
         rna_norm_path = os.path.join(paths['base'], config.RNA_NORMALIZED_FILENAME)
         try:
             rna_normalized.to_csv(rna_norm_path, index=True)
             print(f"Saved normalized RNA FEATURES to {rna_norm_path}")
         except Exception as e: print(f"Warning: Could not save normalized RNA data: {e}")
    # Optionally save final RNA metadata
    # if rna_metadata_final is not None and save_files:
    #     rna_meta_final_path = os.path.join(paths['base'], 'rna_metadata_final_indexed.csv')
    #     rna_metadata_final.to_csv(rna_meta_final_path, index=True)
    #     print(f"Saved final RNA METADATA to {rna_meta_final_path}")

else:
    print("Skipping RNA normalization - indexed RNA features not available.")


# --- 5. Align Samples ---
print("\n--- 5. Aligning Samples ---")
# *** Use the FEATURE dataframes for alignment ***
metabolite_matched, rna_matched = preprocess.align_samples(
    metabolite_features, # Use features output from step 3 (with composite index)
    rna_normalized,      # Use normalized features output from step 4 (with composite index)
    df1_name="Metabolites", df2_name="RNA"
)
# Save matched FEATURE data
if metabolite_matched is not None and rna_matched is not None:
    print(f"Alignment successful. Matched {len(metabolite_matched)} samples.")
    if save_files:
         metab_match_path = os.path.join(paths['base'], config.MATCHED_METABOLITE_FILENAME)
         rna_match_path = os.path.join(paths['base'], config.MATCHED_RNA_FILENAME)
         try:
             metabolite_matched.to_csv(metab_match_path, index=True)
             rna_matched.to_csv(rna_match_path, index=True)
             print(f"Saved matched FEATURE data to {paths['base']}")
         except Exception as e: print(f"Warning: Could not save matched data: {e}")
else:
     print("ERROR: Sample alignment failed. Cannot proceed with downstream analysis.")
     sys.exit(1) # Exit if alignment failed


# --- 6. Variance Filtering ---
print("\n--- 6. Variance Filtering ---")
# Filter the MATCHED feature dataframes
rna_data_filtered = preprocess.apply_variance_filter(
    rna_matched, config.VAR_FILTER_TOP_N_GENES, "RNA"
)
metabolite_data_filtered = preprocess.apply_variance_filter(
    metabolite_matched, config.VAR_FILTER_TOP_N_METABOLITES, "Metabolite"
)
if rna_data_filtered is None or metabolite_data_filtered is None:
     print("ERROR: Variance filtering failed. Cannot proceed.")
     sys.exit(1)


# --- 7. Spearman Correlation ---
print("\n--- 7. Spearman Correlation ---")
# Run on the variance filtered feature dataframes
df_corr_raw = analysis.run_spearman_correlation(rna_data_filtered, metabolite_data_filtered)
# Save raw correlations
if df_corr_raw is not None and save_files:
     corr_raw_path = os.path.join(paths['spearman'], config.CORRELATION_RAW_FILENAME)
     try:
         df_corr_raw.to_csv(corr_raw_path, index=False)
         print(f"Saved raw Spearman results to {corr_raw_path}")
     except Exception as e: print(f"Warning: Could not save raw correlation results: {e}")
elif df_corr_raw is None:
     print("ERROR: Spearman correlation failed. Cannot proceed.")
     sys.exit(1)


# --- 8. P-value Adjustment ---
print("\n--- 8. P-value Adjustment ---")
df_corr_adj = postprocess.adjust_pvalues_bh(df_corr_raw) # Takes raw results df
# Save adjusted correlations
if df_corr_adj is not None and save_files:
     corr_adj_path = os.path.join(paths['spearman'], config.CORRELATION_ADJ_FILENAME)
     try:
         df_corr_adj.to_csv(corr_adj_path, index=False)
         print(f"Saved adjusted Spearman results to {corr_adj_path}")
     except Exception as e: print(f"Warning: Could not save adjusted correlation results: {e}")
elif df_corr_adj is None: # If adjustment itself failed
     print("ERROR: P-value adjustment failed.")
     df_corr_adj = df_corr_raw # Fallback to raw for plotting


# --- 9. Correlation Visualization ---
print("\n--- 9. Correlation Visualization ---")
if df_corr_adj is not None: # Use adjusted (or raw fallback) results
     print(f"Generating correlation plots (saving to {paths['spearman_plots']})...")
     visualize.plot_rho_distribution(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_volcano(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_correlation_clustermap(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_top_correlation_heatmap(df_corr_adj, save_path=paths['spearman_plots'])
     # Scatter plots need the VARIANCE FILTERED feature data
     visualize.plot_top_scatter(df_corr_adj, rna_data_filtered, metabolite_data_filtered, save_path=paths['spearman_plots'])
else:
     print("Skipping correlation visualization as results are not available.")


# --- 10. NMF Analysis & Evaluation ---
print("\n--- 10. NMF Analysis & Evaluation ---")
# Example for multiple k values
k_values_to_test = [3, 4, 5, 6, 7] # Example range
all_nmf_metrics = {} # Dictionary to store metrics for each k

for k_to_run in k_values_to_test:
    print("\n" + "*"*30); print(f" Running NMF for k = {k_to_run}"); print("*"*30)

    # Ensure VARIANCE FILTERED data is available
    if rna_data_filtered is None or metabolite_data_filtered is None:
        print(f"Skipping NMF for k={k_to_run}: Filtered data not available.")
        continue

    nmf_results = analysis.run_nmf_concatenated(
        rna_data_filtered,          # Use variance filtered features
        metabolite_data_filtered,   # Use variance filtered features
        n_components=k_to_run,
        max_iter=config.NMF_DEFAULT_MAX_ITER,
        random_state=config.NMF_DEFAULT_RANDOM_STATE
    )

    if nmf_results["successful"]:
        # Save NMF results
        if save_files:
             try: # Save files
                 h_filename = config.NMF_H_FILENAME_TEMPLATE.format(k=k_to_run)
                 w_rna_filename = config.NMF_W_RNA_FILENAME_TEMPLATE.format(k=k_to_run)
                 w_metab_filename = config.NMF_W_METAB_FILENAME_TEMPLATE.format(k=k_to_run)
                 # Check if DFs exist before saving
                 if nmf_results["H_df"] is not None: nmf_results["H_df"].to_csv(os.path.join(paths['nmf'], h_filename), index=True)
                 if nmf_results["W_rna_df"] is not None: nmf_results["W_rna_df"].to_csv(os.path.join(paths['nmf'], w_rna_filename), index=True)
                 if nmf_results["W_metab_df"] is not None: nmf_results["W_metab_df"].to_csv(os.path.join(paths['nmf'], w_metab_filename), index=True)
                 print(f"Saved NMF results for k={k_to_run} to {paths['nmf']}")
             except Exception as e: print(f"Warning: Could not save NMF results for k={k_to_run}: {e}")

        # Evaluate NMF results
        print(f"\n--- Evaluating NMF k={k_to_run} ---")
        # Evaluation function needs the results components and original (filtered) data for V recalc
        evaluation_metrics = visualize.evaluate_nmf(
            nmf_results["H_df"],
            nmf_results["W_rna_df"],
            nmf_results["W_metab_df"],
            rna_data_filtered=rna_data_filtered, # Pass original filtered data for V recalc
            metabolite_data_filtered=metabolite_data_filtered,
            model=nmf_results["model"],
            k=k_to_run,
            plots_save_path=paths['nmf_plots'] # Pass the specific NMF plots path
        )
        # Store metrics for later comparison
        if evaluation_metrics:
             all_nmf_metrics[k_to_run] = evaluation_metrics
             print(f"Stored evaluation metrics for k={k_to_run}")
    else:
        print(f"NMF analysis failed for k={k_to_run}.")

# --- Post-NMF Evaluation (Example: Plot Reconstruction Error) ---
# *** Initialize lists BEFORE the 'if' check ***
k_vals = []
recon_errors = []
relative_errors = []
plot_k_recon, plot_recon_errors = [], []
plot_k_rel, plot_relative_errors = [], []
plot_relative_errors_aligned = []
plot_k_rel_aligned = []


if all_nmf_metrics: # Check if the dictionary has entries
     print("\n--- NMF Post-Evaluation: Reconstruction Error vs. k ---")
     k_vals = sorted(all_nmf_metrics.keys()) # Now we overwrite if metrics exist

     # Safely extract metrics using .get with default NaN
     recon_errors = [all_nmf_metrics[k].get('Reconstruction Error', np.nan) for k in k_vals]
     relative_errors = [all_nmf_metrics[k].get('Relative Error', np.nan) for k in k_vals]

     # Filter out k values where error calculation failed before plotting
     valid_recon = [(k, e) for k, e in zip(k_vals, recon_errors) if pd.notna(e)]
     valid_rel = [(k, e) for k, e in zip(k_vals, relative_errors) if pd.notna(e)]

     # Unzip the filtered lists if they are not empty
     plot_k_recon, plot_recon_errors = zip(*valid_recon) if valid_recon else ([], [])
     plot_k_rel, plot_relative_errors = zip(*valid_rel) if valid_rel else ([], [])

     # Check if we have any valid reconstruction errors to plot
     if plot_k_recon: # Check if list is not empty
         fig, ax1 = plt.subplots(figsize=(8, 5))
         color = 'tab:red'
         ax1.set_xlabel('Number of Components (k)')
         ax1.set_ylabel('Reconstruction Error (Frobenius)', color=color)
         ax1.plot(plot_k_recon, plot_recon_errors, marker='o', color=color, label='Reconstruction Error') # Use filtered lists
         ax1.tick_params(axis='y', labelcolor=color)
         ax1.grid(True, axis='x', linestyle=':')

         # Only plot relative error if valid values exist
         if plot_k_rel:
             # Align relative errors with the k values plotted for reconstruction error
             plot_relative_errors_aligned = [rel for k, rel in valid_rel if k in plot_k_recon]
             plot_k_rel_aligned = [k for k, rel in valid_rel if k in plot_k_recon]

             if plot_k_rel_aligned: # Check if any aligned relative errors exist
                 ax2 = ax1.twinx()
                 color = 'tab:blue'
                 ax2.set_ylabel('Relative Error', color=color)
                 # Plot only the aligned relative errors
                 ax2.plot(plot_k_rel_aligned, plot_relative_errors_aligned, marker='s', linestyle='--', color=color, label='Relative Error')
                 ax2.tick_params(axis='y', labelcolor=color)

         plt.title('NMF Reconstruction Error vs. Number of Components (k)')
         # Set x-axis ticks explicitly to integer k values if possible
         # Use k_vals used for plotting, which is plot_k_recon here
         ax1.set_xticks(plot_k_recon if plot_k_recon else k_vals) # Use the k's actually plotted if possible
         fig.tight_layout()
         # Save this plot too
         # Use the correct plots path from the 'paths' dictionary
         plot_filename_k_eval = os.path.join(paths['nmf_plots'], "nmf_k_evaluation_recon_error.png")
         try:
             plt.savefig(plot_filename_k_eval, dpi=150, bbox_inches='tight')
             print(f"Saved k evaluation plot to: {plot_filename_k_eval}")
         except Exception as e_p: print(f"Error saving k eval plot: {e_p}")
         plt.show()
         plt.close()
     else:
         print("Could not plot reconstruction errors (no valid error values found in metrics dictionary).")

else:
    print("\nNo NMF evaluation metrics collected (NMF might have failed for all k or dictionary empty).") # Updated message


print("\n--- FULL ANALYSIS PIPELINE COMPLETE ---")