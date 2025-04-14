# Example Script using the MeERCAT package (Local File Version)

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
    # Assumes the script is in a subdirectory (like 'examples')
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
    # *** UPDATED IMPORT STATEMENTS ***
    import meercat
    from meercat import config
    from meercat import utils
    from meercat import load_data
    from meercat import preprocess
    from meercat import analysis
    from meercat import postprocess
    from meercat import visualize
except ImportError as e:
     print(f"ERROR: Could not import the 'meercat' package.") # Updated name
     print(f"Ensure the package is installed (`pip install .` in the root directory)")
     print(f"or the path is correctly added to sys.path if running from source.")
     print(f"Import error details: {e}")
     sys.exit(1) # Exit if package not found

# ==================================================
# --- User Configuration ---
# ==================================================
# --- Define Base Path for Input Data and Output Results ---
# --- *** USER MUST SET THIS PATH *** ---
# This is the main folder where 'input_data' lives and 'output_results' will be created.
# Example: '/Users/yourname/Documents/MyMultiOmicsProject/' or 'C:/Users/yourname/MyProject/'
# Use '.' for the current directory where the script is run.
BASE_PROJECT_PATH = '.'

# --- Define Subdirectory Names (can use defaults from config or customize) ---
INPUT_DATA_SUBDIR = 'input_data'   # Expects subfolders 'metabolites', 'rna' inside this
OUTPUT_RESULTS_SUBDIR = 'output_results' # Where all analysis results will be saved
# ==================================================


# --- 1. Setup Paths ---
print("--- 1. Setting up Paths ---")
# Construct full paths based on user config
input_data_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, INPUT_DATA_SUBDIR))
output_results_path = os.path.abspath(os.path.join(BASE_PROJECT_PATH, OUTPUT_RESULTS_SUBDIR))

# Use the setup_paths utility from the package to create output structure
# It uses defaults from config.py for subfolder names unless overridden
# Pass the desired BASE output path
paths = utils.setup_paths(
    output_results_path, # Base path for *outputs*
    '',                  # Project folder name within output (use empty if base IS output folder)
    config.DEFAULT_SPEARMAN_SUBFOLDER,
    config.DEFAULT_NMF_SUBFOLDER,
    config.DEFAULT_PLOTS_SUBFOLDER # This plots subfolder is relative inside spearman/nmf
)
# Add specific input paths for clarity
paths['input'] = input_data_path
paths['input_metabolites'] = os.path.join(input_data_path, 'metabolites')
paths['input_rna'] = os.path.join(input_data_path, 'rna')
paths['input_metadata'] = os.path.join(input_data_path, 'metadata.csv') # Assumes metadata is here

print(f"\nInput Data expected in: {paths['input']}")
print(f"Output Results will be saved under: {output_results_path}")
print(f"  Spearman Results: {paths['spearman']}")
print(f"  NMF Results: {paths['nmf']}")
print(f"  Spearman Plots: {paths['spearman_plots']}")
print(f"  NMF Eval Plots: {paths['nmf_plots']}")

# Flag for saving intermediate/final files
save_files = True # Set to False to only run analysis in memory


# --- 2. Load Data ---
print("\n--- 2. Loading Data ---")
# --- Load Metadata ---
metadata_df = load_data.load_metadata(paths['input_metadata'])
if metadata_df is None:
     print("WARNING: Metadata failed to load. Some processing steps might be affected.")
     # Decide whether to exit or continue without metadata if possible
     # sys.exit(1)

# --- Load Metabolites ---
# Check if metabolite input directory exists
if not os.path.isdir(paths['input_metabolites']):
     print(f"ERROR: Metabolite input directory not found: {paths['input_metabolites']}")
     loaded_metabolite_data = {} # Set to empty
else:
     metabolite_file_dict = {f: os.path.join(paths['input_metabolites'], f)
                              for f in os.listdir(paths['input_metabolites']) if f.lower().endswith('.csv')}
     loaded_metabolite_data = load_data.load_metabolite_files(metabolite_file_dict)
     if not loaded_metabolite_data:
          print("ERROR: No metabolite files loaded. Check directory contents and file formats.")
          # sys.exit(1)

# --- Load RNA ---
if not os.path.isdir(paths['input_rna']):
     print(f"ERROR: RNA input directory not found: {paths['input_rna']}")
     loaded_rna_data = {} # Set to empty
else:
     rna_file_dict = {f: os.path.join(paths['input_rna'], f)
                       for f in os.listdir(paths['input_rna']) if f.lower().endswith('.csv')} # Adjust if needed
     # *** IMPORTANT: Set rows_are_genes correctly for your data format ***
     loaded_rna_data = load_data.load_rna_files(rna_file_dict, rows_are_genes=True)
     if not loaded_rna_data:
          print("ERROR: No RNA files loaded. Check directory contents and file formats.")
          # sys.exit(1)


# --- 3. Preprocess Metabolites ---
print("\n--- 3. Preprocessing Metabolites ---")
metabolite_cleaned_indexed = None
metabolite_combined_raw = preprocess.combine_dataframes(loaded_metabolite_data, axis=0)
if metabolite_combined_raw is not None:
    metabolite_cleaned_indexed = preprocess.clean_metabolite_data(metabolite_combined_raw)
    # Save cleaned data (optional intermediate step)
    if metabolite_cleaned_indexed is not None and save_files:
        metab_clean_path = os.path.join(paths['base'], config.METABOLITE_CLEANED_FILENAME)
        try:
            metabolite_cleaned_indexed.to_csv(metab_clean_path, index=True)
            print(f"Saved cleaned metabolite data to {metab_clean_path}")
        except Exception as e: print(f"Warning: Could not save cleaned metabolite data: {e}")
else: print("Skipping metabolite preprocessing - no data loaded.")


# --- 4. Preprocess RNA ---
print("\n--- 4. Preprocessing RNA ---")
rna_normalized = None
processed_rna_dict = {}
if loaded_rna_data: # Check if dict has content
    for name, df_raw in loaded_rna_data.items():
        # *** IMPORTANT: Set rows_are_genes correctly for your data format ***
        processed_df = preprocess.process_rna_dataframe(df_raw, name, rows_are_genes=True)
        if processed_df is not None:
            processed_rna_dict[name] = processed_df
rna_combined_processed = preprocess.combine_dataframes(processed_rna_dict, axis=0, join='outer')
if rna_combined_processed is not None:
    rna_normalized = preprocess.normalize_rna_data(rna_combined_processed)
    # Save normalized data (optional intermediate step)
    if rna_normalized is not None and save_files:
         rna_norm_path = os.path.join(paths['base'], config.RNA_NORMALIZED_FILENAME)
         try:
             rna_normalized.to_csv(rna_norm_path, index=True)
             print(f"Saved normalized RNA data to {rna_norm_path}")
         except Exception as e: print(f"Warning: Could not save normalized RNA data: {e}")
else: print("Skipping RNA preprocessing - no data loaded.")


# --- 5. Align Samples ---
print("\n--- 5. Aligning Samples ---")
metabolite_matched, rna_matched = preprocess.align_samples(
    metabolite_cleaned_indexed, # Use output from step 3
    rna_normalized,             # Use output from step 4
    df1_name="Metabolites", df2_name="RNA"
)
# Save matched data
if metabolite_matched is not None and rna_matched is not None:
    print(f"Alignment successful. Matched {len(metabolite_matched)} samples.")
    if save_files:
         metab_match_path = os.path.join(paths['base'], config.MATCHED_METABOLITE_FILENAME)
         rna_match_path = os.path.join(paths['base'], config.MATCHED_RNA_FILENAME)
         try:
             metabolite_matched.to_csv(metab_match_path, index=True)
             rna_matched.to_csv(rna_match_path, index=True)
             print(f"Saved matched data to {paths['base']}")
         except Exception as e: print(f"Warning: Could not save matched data: {e}")
else:
     print("ERROR: Sample alignment failed. Cannot proceed with downstream analysis.")
     sys.exit(1) # Exit if alignment failed


# --- 6. Variance Filtering ---
print("\n--- 6. Variance Filtering ---")
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
df_corr_adj = postprocess.adjust_pvalues_bh(df_corr_raw)
# Save adjusted correlations
if df_corr_adj is not None and save_files:
     corr_adj_path = os.path.join(paths['spearman'], config.CORRELATION_ADJ_FILENAME)
     try:
         df_corr_adj.to_csv(corr_adj_path, index=False)
         print(f"Saved adjusted Spearman results to {corr_adj_path}")
     except Exception as e: print(f"Warning: Could not save adjusted correlation results: {e}")
else:
     print("ERROR: P-value adjustment failed. Correlation visualization might use raw p-values or fail.")
     # df_corr_adj = df_corr_raw # Fallback to raw if needed


# --- 9. Correlation Visualization ---
print("\n--- 9. Correlation Visualization ---")
if df_corr_adj is not None:
     print(f"Generating correlation plots (saving to {paths['spearman_plots']})...")
     visualize.plot_rho_distribution(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_volcano(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_correlation_clustermap(df_corr_adj, save_path=paths['spearman_plots'])
     visualize.plot_top_correlation_heatmap(df_corr_adj, save_path=paths['spearman_plots'])
     # Scatter plots need the filtered data passed explicitly
     visualize.plot_top_scatter(df_corr_adj, rna_data_filtered, metabolite_data_filtered, save_path=paths['spearman_plots'])
else:
     print("Skipping correlation visualization as adjusted results are not available.")


# --- 10. NMF Analysis & Evaluation ---
print("\n--- 10. NMF Analysis & Evaluation ---")
# Example for multiple k values
k_values_to_test = [3, 4, 5, 6, 7] # Example range
all_nmf_metrics = {} # Dictionary to store metrics for each k

for k_to_run in k_values_to_test:
    print("\n" + "*"*30)
    print(f" Running NMF for k = {k_to_run}")
    print("*"*30)

    # Make sure filtered data is available for NMF run
    if rna_data_filtered is None or metabolite_data_filtered is None:
        print(f"Skipping NMF for k={k_to_run}: Filtered data not available.")
        continue # Skip to next k

    nmf_results = analysis.run_nmf_concatenated(
        rna_data_filtered,
        metabolite_data_filtered,
        n_components=k_to_run, # Use loop variable
        max_iter=config.NMF_DEFAULT_MAX_ITER,
        random_state=config.NMF_DEFAULT_RANDOM_STATE
    )

    if nmf_results["successful"]:
        # Save NMF results
        if save_files:
             try:
                 h_filename = config.NMF_H_FILENAME_TEMPLATE.format(k=k_to_run)
                 w_rna_filename = config.NMF_W_RNA_FILENAME_TEMPLATE.format(k=k_to_run)
                 w_metab_filename = config.NMF_W_METAB_FILENAME_TEMPLATE.format(k=k_to_run)
                 nmf_results["H_df"].to_csv(os.path.join(paths['nmf'], h_filename), index=True)
                 nmf_results["W_rna_df"].to_csv(os.path.join(paths['nmf'], w_rna_filename), index=True)
                 nmf_results["W_metab_df"].to_csv(os.path.join(paths['nmf'], w_metab_filename), index=True)
                 print(f"Saved NMF results for k={k_to_run} to {paths['nmf']}")
             except Exception as e: print(f"Warning: Could not save NMF results for k={k_to_run}: {e}")

        # Evaluate NMF results
        print(f"\n--- Evaluating NMF k={k_to_run} ---")
        # Evaluation function needs the results components and original (filtered) data for V recalc
        evaluation_metrics = visualize.evaluate_nmf(
            nmf_results["H_df"],
            nmf_results["W_rna_df"],
            nmf_results["W_metab_df"],
            # Pass original filtered data FOR V RECALCULATION inside evaluate_nmf
            # This avoids keeping the large V_combined_imputed in memory for all k
            rna_data_filtered=rna_data_filtered, # Pass original filtered data
            metabolite_data_filtered=metabolite_data_filtered, # Pass original filtered data
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
if all_nmf_metrics:
     print("\n--- NMF Post-Evaluation: Reconstruction Error vs. k ---")
     k_vals = sorted(all_nmf_metrics.keys())
     recon_errors = [all_nmf_metrics[k].get('Reconstruction Error', np.nan) for k in k_vals]
     relative_errors = [all_nmf_metrics[k].get('Relative Error', np.nan) for k in k_vals]

     # Check if we have errors to plot
     if not all(np.isnan(e) for e in recon_errors):
         fig, ax1 = plt.subplots(figsize=(8, 5))
         color = 'tab:red'
         ax1.set_xlabel('Number of Components (k)')
         ax1.set_ylabel('Reconstruction Error (Frobenius)', color=color)
         ax1.plot(k_vals, recon_errors, marker='o', color=color, label='Reconstruction Error')
         ax1.tick_params(axis='y', labelcolor=color)
         ax1.grid(True, axis='x', linestyle=':')

         # Only plot relative error if calculated
         if not all(np.isnan(e) for e in relative_errors):
             ax2 = ax1.twinx()
             color = 'tab:blue'
             ax2.set_ylabel('Relative Error', color=color)
             ax2.plot(k_vals, relative_errors, marker='s', linestyle='--', color=color, label='Relative Error')
             ax2.tick_params(axis='y', labelcolor=color)

         plt.title('NMF Reconstruction Error vs. Number of Components (k)')
         fig.tight_layout()
         # Save this plot too
         plot_filename_k_eval = os.path.join(paths['nmf_plots'], "nmf_k_evaluation_recon_error.png")
         try: plt.savefig(plot_filename_k_eval, dpi=150, bbox_inches='tight'); print(f"Saved k evaluation plot to: {plot_filename_k_eval}")
         except Exception as e_p: print(f"Error saving k eval plot: {e_p}")
         plt.show()
         plt.close()
     else:
         print("Could not plot reconstruction errors (values might be NaN).")

else:
    print("\nNo NMF evaluation metrics collected (NMF might have failed for all k).")


print("\n--- FULL ANALYSIS PIPELINE COMPLETE ---")