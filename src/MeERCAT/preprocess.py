# mercat_analyzer/meercat/preprocess.py

import pandas as pd



## Helper Function 1: Get data dictionary for metabolites
def get_metabolite_data_dict(all_loaded_data):
    """
    Helper function to extract the metabolite data dictionary
    from the main loaded data dictionary.

    Args:
        all_loaded_data: The dictionary returned by load_all_data.

    Returns:
        A dictionary containing the metabolite data
        ({exp_id: DataFrame, ...}), or an empty dictionary {}
        if 'metabolite_data' key is not found in the input.
    """
    # Use .get() to safely retrieve the dictionary.
    # If 'metabolite_data' key doesn't exist, it returns the default value ({}).
    metabolite_dict = all_loaded_data.get("metabolite_data", {})
    return metabolite_dict



## Helper Function 2: Get data dictionary for RNA
def get_rna_data_dict(all_loaded_data):
    """
    Helper function to extract the RNA data dictionary
    from the main loaded data dictionary.

    Args:
        all_loaded_data: The dictionary returned by load_all_data.

    Returns:
        A dictionary containing the RNA data
        ({exp_id: DataFrame, ...}), or an empty dictionary {}
        if 'rna_data' key is not found in the input.
    """
    # Use .get() to safely retrieve the dictionary.
    rna_dict = all_loaded_data.get("rna_data", {})
    return rna_dict



# Helper function 3: Metabolite Merge Function  
def merge_metabolite_data(all_loaded_data) -> pd.DataFrame:
    """
    Extracts metabolite data using a helper and merges the DataFrames
    from the resulting dictionary into a single DataFrame.

    Args:
        all_loaded_data:
            The main dictionary returned by load_all_data, containing potentially
            'metabolite_data' as a key.

    Returns:
        pd.DataFrame: A single DataFrame containing all rows from the extracted
                      metabolite DataFrames concatenated vertically. Returns an empty
                      DataFrame if no metabolite data was found or merged.
    """
    print("Merging Metabolite Data...")
    metabolite_data_dict = get_metabolite_data_dict(all_loaded_data)
    list_of_metabolite_dfs = list(metabolite_data_dict.values())
    # Concatenate the list of DataFrames vertically
    try:
        merged_metabolites = pd.concat(list_of_metabolite_dfs, ignore_index=True)
        print(f"  Merged {len(list_of_metabolite_dfs)} metabolite DataFrames.")
        print(f"  Dataframe shape: {merged_metabolites.shape}")
        print("----------------------------\n\n")
        return merged_metabolites
    except Exception as e:
        print(f"Merge Function: Error during pd.concat: {e}")
        print("Merge Function: Returning empty DataFrame due to concatenation error.")
        print("----------------------------\n\n")
        return pd.DataFrame()


#Helper function 4: RNA Merge Function
def merge_rna_data(all_loaded_data) -> pd.DataFrame:
    """
    Extracts RNA data using a helper and merges the DataFrames
    from the resulting dictionary into a single DataFrame.Then it transposes the DataFrame
    to have samples as rows and genes as columns.

    Args:
        all_loaded_data:
            The main dictionary returned by load_all_data, containing potentially
            'rna_data' as a key.

    Returns:
        pd.DataFrame: A single transformed DataFrame containing all rows from the extracted
                      RNA DataFrames concatenated vertically. Returns an empty
                      DataFrame if no RNA data was found or merged.
    """
    print("Merging RNA Data...")
    rna_data_dict = get_rna_data_dict(all_loaded_data)
    list_of_rna_dfs = list(rna_data_dict.values())
    # Concatenate the list of DataFrames horizontally
    try:
        # Use axis=1 for horizontal merge.
        # Pandas will align based on the index (rows).
        print(f"  Merged {len(list_of_rna_dfs)} RNA DataFrames.")
        merged_rna = pd.concat(list_of_rna_dfs, axis=1)
        print(f"  Dataframe shape befor for transpose: {merged_rna.shape}")
        merged_rna= merged_rna.T
        print(f"  Dataframe shape after transpose: {merged_rna.shape}")
        print("----------------------------\n\n")
    except Exception as e:
        print(f"Merge Function: Error during horizontal pd.concat: {e}")
        print("Merge Function: Returning empty DataFrame due to concatenation error.")
        print("----------------------------\n\n")
        return pd.DataFrame()























def add_sample_id_to_rna(
    merged_rna_transposed: pd.DataFrame,
    rna_metadata: pd.DataFrame,
    sample_id_cols: list = None
    ) -> pd.DataFrame:
    """
    Adds a 'Sample_ID' column to the merged (transposed) RNA DataFrame.

    The new ID is created by combining specified columns from the RNA metadata,
    aligned by the original RNA sample ID (index of merged_rna_transposed).

    Args:
        merged_rna_transposed (pd.DataFrame): The RNA DataFrame with samples
                                             as index and features as columns.
        rna_metadata (pd.DataFrame): The RNA metadata DataFrame, indexed by
                                     original RNA sample ID.
        sample_id_cols (list, optional): List of column names from
                                            rna_metadata to combine for the ID.
                                            Defaults to config.RNA_EXTERNAL_METADATA_COLS.

    Returns:
        pd.DataFrame: The RNA DataFrame with the added 'Sample_ID' column,
                      or the original DataFrame if inputs are invalid or
                      an error occurs.
    """
    print("\n--- Adding Sample ID to RNA Data ---")

    # Determine columns to use for the sample ID
    if sample_id_cols is None:
        try:
            # Default to the list defined in config.py
            sample_id_cols = config.RNA_EXTERNAL_METADATA_COLS
            if not sample_id_cols: # Check if config list itself is empty
                 raise ValueError("Config list RNA_EXTERNAL_METADATA_COLS is empty.")
            print(f"  Using default sample ID columns from config: {sample_id_cols}")
        except (AttributeError, NameError, TypeError, ValueError) as e:
             print(f"  * Error: Could not get valid default sample ID columns from config: {e}")
             print("  Please provide 'sample_id_cols' list argument.")
             return merged_rna_transposed

    if not isinstance(sample_id_cols, list) or not sample_id_cols:
         print("  * Error: 'sample_id_cols' must be a non-empty list of column names.")
         return merged_rna_transposed

    # --- Prepare for Alignment ---
    # Work on copies to avoid modifying original inputs unexpectedly
    rna_data = merged_rna_transposed.copy()
    metadata = rna_metadata.copy()

    # Clean the index of both DataFrames for reliable matching
    # It's crucial that rna_metadata's original index WAS the sample ID
    print("  Cleaning indices for alignment...")
    try:
        original_rna_data_index_name = rna_data.index.name # Store original name if exists
        rna_data.index = clean_index_simple(rna_data.index)

        original_metadata_index_name = metadata.index.name # Store original name if exists
        metadata.index = clean_index_simple(metadata.index)
        print("  Indices cleaned.")
    except Exception as e:
        print(f"  * Error cleaning indices: {e}")
        return merged_rna_transposed # Return original on error

    # Check if required columns exist in RNA metadata *after* cleaning index
    missing_cols = [col for col in sample_id_cols if col not in metadata.columns]
    if missing_cols:
        print(f"  * Error: Required columns for composite ID missing in RNA metadata: {missing_cols}")
        return merged_rna_transposed

    # --- Align and Create ID ---
    try:
        print("  Aligning metadata and creating composite IDs...")
        # Reindex metadata to match the order and presence of indices in rna_data
        # This effectively performs a left-join lookup based on the cleaned index
        aligned_metadata = metadata.loc[rna_data.index, sample_id_cols]

        # Create the composite ID string vectorially
        # Replace potential NaN/None values with 'NA' string before concatenation
        composite_id_series = aligned_metadata[sample_id_cols[0]].fillna('NA').astype(str)
        for col in sample_id_cols[1:]:
            composite_id_series += '_' + aligned_metadata[col].fillna('NA').astype(str)

        # --- Add Composite ID as a Column ---
        # Assign the generated Series to the new column.
        # Since aligned_metadata was indexed by rna_data.index, the series aligns correctly.
        rna_data['Sample_ID'] = composite_id_series.values
        print("  Successfully added 'Sample_ID' column.")

        # Restore original index name if desired (optional)
        # rna_data.index.name = original_rna_data_index_name

        return rna_data # Return the modified DataFrame

    except KeyError as e:
         print(f"  * Error: A key error occurred during alignment (likely mismatch between RNA data index and metadata index): {e}")
         print(f"  Check if the indices contain the same sample IDs after cleaning.")
         return merged_rna_transposed
    except Exception as e:
        print(f"  * Error during alignment or composite ID creation: {e}")
        return merged_rna_transposed







# import pandas as pd
# import numpy as np
# import re
# # Import utils and config functions/variables
# from .utils import identify_metabolite_columns, clean_feature_names, handle_duplicate_features
# from .config import (
#     METADATA_SAMPLE_COL, SAMPLE_FILTER_EXCLUDE_PATTERNS,
#     # METADATA_sample_id_cols, # No longer needed for index creation here
#     RNA_FILTER_MIN_COUNTS,
#     RNA_FILTER_MIN_SAMPLES_FRAC, VAR_FILTER_TOP_N_GENES,
#     VAR_FILTER_TOP_N_METABOLITES, METABOLITE_ID_PREFIXES,
#     METABOLITE_ID_SUFFIX_REGEX #, RNA_METADATA_COLS_EXPECTED # Not used directly here now
# )
# # Import clean_index from load_data
# try:
#     from .load_data import clean_index
# except ImportError:
#     print("ERROR in preprocess.py: Could not import clean_index from .load_data.")
#     def clean_index(index): print("ERROR: clean_index function missing!"); return index


# # --- combine_dataframes function ---
# def combine_dataframes(df_dict, axis=0, join='outer'):
#     """Combines a dictionary of DataFrames."""
#     if not isinstance(df_dict, dict) or not df_dict: print("Warning: Input not valid dict."); return None
#     print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
#     try:
#         valid_dfs = [df for df in df_dict.values() if isinstance(df, pd.DataFrame) and not df.empty]
#         if not valid_dfs: print("Warning: No valid DataFrames found."); return None
#         # Preserve index when combining rows (axis=0)
#         combined_df = pd.concat(valid_dfs, axis=axis, join=join, ignore_index=False) # ignore_index=False
#         print(f"Combined shape: {combined_df.shape}")
#         # Use INFO for duplicate index warning as it's expected when combining raw files
#         if axis == 0 and combined_df.index.has_duplicates: print(f"INFO: Input indices had duplicates (expected for axis=0 combine).")
#         return combined_df
#     except Exception as e: print(f"Error during concatenation: {e}"); return None


# def clean_metabolite_data(combined_df):
#     """Cleans combined metabolite data, averages duplicates, separates metadata/features,
#        filters samples, and sets index using the SHORT version of 'original_rna_sample_id'.
#     """
#     print("\n--- Cleaning Combined Metabolite Data ---")
#     if combined_df is None or combined_df.empty: return None, None
#     df_cleaned = combined_df.copy(); metadata_df_final = None; df_features_only = None; index_set_successfully = False
#     RAW_ID_COL_NAME = 'original_rna_sample_id'
#     FINAL_INDEX_NAME = 'original_rna_sample_id' # Keep the final desired name

#     try:
#         if RAW_ID_COL_NAME not in df_cleaned.columns: raise ValueError(f"Column '{RAW_ID_COL_NAME}' not found.")

#         # --- Create SHORT sample ID (without prefix) FIRST ---
#         print(f"Extracting short sample ID from '{RAW_ID_COL_NAME}'...")
#         short_id_series = df_cleaned[RAW_ID_COL_NAME].astype(str).str.extract(r'(?:mc\d+_)?(.*)', expand=False)
#         short_id_series = clean_index(short_id_series) # Clean the extracted short ID
#         SHORT_ID_COL_TEMP = 'short_sample_id_temp' # Temporary column name
#         df_cleaned[SHORT_ID_COL_TEMP] = short_id_series

#         # --- Handle Duplicates based on the SHORT ID ---
#         print(f"Checking for duplicates in generated short sample ID...")
#         duplicates_present = df_cleaned[SHORT_ID_COL_TEMP].duplicated().any()
#         index_is_unique = False

#         if duplicates_present:
#             print(f"  Duplicates found in short ID. Averaging numeric features...")
#             feature_cols_initial, non_feature_cols_initial = identify_metabolite_columns(df_cleaned.columns.drop(SHORT_ID_COL_TEMP), METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX)
#             numeric_feature_cols = df_cleaned[feature_cols_initial].select_dtypes(include=np.number).columns.tolist()
#             metadata_cols_to_keep = [col for col in non_feature_cols_initial if col in df_cleaned.columns] # Keep original metadata cols
#             # Group by the SHORT ID now
#             grouping_cols = [SHORT_ID_COL_TEMP] + metadata_cols_to_keep
#             agg_dict = {**{col: 'first' for col in metadata_cols_to_keep}, **{col: 'mean' for col in numeric_feature_cols}}
#             # Group by the short ID, keep first metadata, average numeric
#             df_grouped = df_cleaned.groupby(SHORT_ID_COL_TEMP).agg(agg_dict)
#             if df_grouped.index.duplicated().any(): raise ValueError("Duplicates remain after averaging!")
#             df_cleaned = df_grouped.reset_index() # SHORT_ID_COL_TEMP becomes a column again
#             print(f"  Shape after averaging duplicates based on short ID: {df_cleaned.shape}")
#             index_is_unique = True
#         else:
#             print(f"  Values based on short sample ID are unique.")
#             index_is_unique = True # Already unique based on short ID

#         # --- Continue with cleaning on the (potentially averaged) df_cleaned ---
#         # (Identify features, Convert numeric, Drop NaNs - as before)
#         all_cols = df_cleaned.columns.tolist()
#         metabolite_feature_cols, non_metabolite_cols = identify_metabolite_columns(all_cols, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX)
#         potential_metadata_cols = list(set([col for col in non_metabolite_cols if col != SHORT_ID_COL_TEMP] + [METADATA_SAMPLE_COL])) # Use temp col name
#         metadata_cols_present = [col for col in potential_metadata_cols if col in df_cleaned.columns]
#         feature_cols_present = [col for col in metabolite_feature_cols if col in df_cleaned.columns]
#         if feature_cols_present:
#             df_features = df_cleaned[feature_cols_present].copy()
#             for col in feature_cols_present:
#                  if not pd.api.types.is_numeric_dtype(df_features[col]): df_features[col] = pd.to_numeric(df_features[col], errors='coerce')
#             numeric_feature_cols = df_features.select_dtypes(include=np.number).columns.tolist()
#             df_numeric_features = df_features[numeric_feature_cols].copy()
#             nan_feature_cols = [col for col in numeric_feature_cols if df_numeric_features[col].isnull().all()]
#             if nan_feature_cols: df_numeric_features.drop(columns=nan_feature_cols, inplace=True); print(f"Dropped {len(nan_feature_cols)} all-NaN features.")
#         else: df_numeric_features = pd.DataFrame(index=df_cleaned.index)

#         # Prepare for sample filtering (re-attach relevant metadata)
#         cols_for_meta_df = [SHORT_ID_COL_TEMP] + metadata_cols_present # Use temp col name
#         metadata_df = df_cleaned[cols_for_meta_df].copy()
#         df_numeric_features.index = metadata_df.index # Align indices before merge
#         df_to_filter = pd.merge(metadata_df, df_numeric_features, left_index=True, right_index=True, how='left')

#         # Filter out QC/Blank samples
#         print("Filtering out QC/Blank samples...")
#         df_filtered = df_to_filter
#         if METADATA_SAMPLE_COL in df_filtered.columns:
#             exclude_regex = '|'.join(SAMPLE_FILTER_EXCLUDE_PATTERNS)
#             keep_mask = ~df_filtered[METADATA_SAMPLE_COL].astype(str).str.lower().str.contains(exclude_regex, na=False, regex=True)
#             rows_before = len(df_filtered)
#             df_filtered = df_filtered[keep_mask].copy()
#             print(f"Filtered {rows_before - len(df_filtered)} QC/Blank rows. Shape now: {df_filtered.shape}")
#         if df_filtered.empty: print("Warning: DataFrame empty after filtering."); return pd.DataFrame(), pd.DataFrame()

#         # Separate final metadata and features AGAIN
#         final_metadata_cols = [col for col in metadata_cols_present if col in df_filtered.columns]
#         metadata_df_final = df_filtered[[SHORT_ID_COL_TEMP] + final_metadata_cols].copy() # Keep temp ID col
#         final_feature_cols = [col for col in df_filtered.columns if col not in final_metadata_cols and col != SHORT_ID_COL_TEMP]
#         df_features_only = df_filtered[final_feature_cols].copy()
#         df_features_only[SHORT_ID_COL_TEMP] = df_filtered[SHORT_ID_COL_TEMP] # Add temp ID col
#         print(f"Separated final Metadata ({metadata_df_final.shape}) and Features ({df_features_only.shape})")

#         # Set index using the SHORT sample ID IF IT WAS UNIQUE
#         if index_is_unique:
#              print(f"Setting index to the derived short sample ID ('{FINAL_INDEX_NAME}')...")
#              try:
#                   metadata_df_final = metadata_df_final.set_index(SHORT_ID_COL_TEMP)
#                   df_features_only = df_features_only.set_index(SHORT_ID_COL_TEMP)
#                   metadata_df_final.index.name = FINAL_INDEX_NAME # Use the desired final name
#                   df_features_only.index.name = FINAL_INDEX_NAME # Use the desired final name
#                   index_set_successfully = True; print("Index set successfully.")
#              except Exception as e_index: print(f"ERROR setting index: {e_index}"); index_set_successfully = False
#         else: print(f"Skipping setting index due to duplicate short IDs generated earlier."); index_set_successfully = False

#     except Exception as e_clean: print(f"ERROR during metabolite cleaning: {e_clean}"); return None, None

#     print("Metabolite cleaning finished.")
#     if index_set_successfully: return df_features_only, metadata_df_final
#     else: print("Warning: Returning metabolite data without index set."); return df_features_only, metadata_df_final



# # --- process_rna_dataframe function ---
# def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
#     """Orients, cleans features, returns RNA features indexed by original sample ID."""
#     print(f"\nProcessing RNA DataFrame from '{filename}'...")
#     if df_raw is None or df_raw.empty: return None
#     features_df = None
#     try:
#         if rows_are_genes:
#             df_reset = df_raw.reset_index(); gene_id_col = df_reset.columns[0]
#             sample_pattern = r'^\d+_\d{2,}_\d+[a-zA-Z]+$'
#             sample_cols = [col for col in df_reset.columns if col != gene_id_col and isinstance(col, str) and re.match(sample_pattern, col)]
#             annotation_cols = [col for col in df_reset.columns if col != gene_id_col and col not in sample_cols]
#             if not sample_cols: sample_cols = [col for col in df_reset.columns if col != gene_id_col]; print(f"  Warning: Using fallback for RNA sample columns in {filename}.")
#             # if annotation_cols: print(f"  Identified {len(annotation_cols)} annotation columns: {annotation_cols}") # Less verbose
#             df_numeric_part = df_reset.set_index(gene_id_col)[sample_cols].apply(pd.to_numeric, errors='coerce')
#             features_df = df_numeric_part.T.copy()
#         else: features_df = df_raw.select_dtypes(include=np.number).apply(pd.to_numeric, errors='coerce').copy()
#         print(f"  Oriented features shape (Samples x Genes): {features_df.shape}")
#         features_df.index = clean_index(features_df.index) # Uses imported clean_index
#         features_df.index.name = 'original_rna_sample_id' # Set index name
#         features_df.columns = clean_index(features_df.columns); features_df.columns = clean_feature_names(features_df.columns); features_df = handle_duplicate_features(features_df)
#         nan_count = features_df.isna().sum().sum()
#         if nan_count > 0: features_df.fillna(0, inplace=True)
#         print(f"Finished basic processing for RNA file {filename}.")
#         return features_df
#     except Exception as e: print(f"  Error during basic RNA processing for '{filename}': {e}"); return None


# # --- normalize_rna_data function ---
# def normalize_rna_data(rna_features_orig_index):
#     """Filters low counts, normalizes (CPM+Log2). Input indexed by original_rna_sample_id."""
#     print("\n--- Filtering and Normalizing Combined RNA Data ---")
#     if rna_features_orig_index is None or rna_features_orig_index.empty: return None
#     if rna_features_orig_index.index.name != 'original_rna_sample_id': print("ERROR: Input must be indexed by 'original_rna_sample_id'."); return None
#     rna_counts = rna_features_orig_index.select_dtypes(include=np.number).fillna(0).astype(np.float32)
#     if rna_counts.empty: print("Error: No numeric gene columns."); return None
#     print(f"Input shape: {rna_counts.shape}")
#     min_samples = max(2, int(len(rna_counts) * RNA_FILTER_MIN_SAMPLES_FRAC))
#     print(f"Filtering: Keeping genes >= {RNA_FILTER_MIN_COUNTS} counts in >= {min_samples} samples.")
#     try: genes_to_keep_mask = (rna_counts.values >= RNA_FILTER_MIN_COUNTS).sum(axis=0) >= min_samples; rna_filtered = rna_counts.loc[:, genes_to_keep_mask]
#     except Exception as e: print(f"Error during gene filtering mask: {e}"); return None
#     if rna_filtered.empty: print("Error: Filtering removed all genes."); return None
#     print(f"Removed {rna_counts.shape[1] - rna_filtered.shape[1]} genes. Shape after filtering: {rna_filtered.shape}")
#     print("Normalizing (CPM + Log2)...")
#     try:
#         library_sizes_col = rna_filtered.values.sum(axis=1)[:, np.newaxis] + 1e-9
#         cpm_array = (rna_filtered.values / library_sizes_col) * 1e6
#         log2_cpm_array = np.log2(np.maximum(cpm_array, 0) + 1)
#         rna_normalized_log2 = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
#         print("Normalization complete. Final shape:", rna_normalized_log2.shape)
#         return rna_normalized_log2
#     except Exception as e: print(f"Error during normalization: {e}"); return None


# # --- align_samples function (Expects original_rna_sample_id index) ---
# def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
#     """Aligns two DataFrames based on common indices (expects 'original_rna_sample_id')."""
#     print(f"\n--- Aligning {df1_name} and {df2_name} to Common Samples ---")
#     if df1 is None or df2 is None or df1.empty or df2.empty: print("Input DataFrame missing/empty."); return None, None
#     if df1.index.name != 'original_rna_sample_id' or df2.index.name != 'original_rna_sample_id': print("Error: Inputs must be indexed by 'original_rna_sample_id'."); return None, None
#     print(f"Input {df1_name} shape: {df1.shape}; Input {df2_name} shape: {df2.shape}")
#     common_samples = df1.index.intersection(df2.index)
#     print(f"Found {len(common_samples)} common samples based on original_rna_sample_id index.")
#     if len(common_samples) == 0: print("CRITICAL WARNING: No common samples found!"); return None, None
#     try:
#         df1_matched = df1.loc[common_samples].copy()
#         df2_aligned = df2.loc[common_samples].reindex(df1_matched.index).copy()
#         if not df1_matched.index.equals(df2_aligned.index): raise ValueError("Final indices mismatch!")
#         print(f"Aligned shapes: {df1_name} {df1_matched.shape}, {df2_name} {df2_aligned.shape}. Alignment successful.")
#         return df1_matched, df2_aligned # Returns indexed by original_rna_sample_id
#     except Exception as e: print(f"Error during sample alignment: {e}"); return None, None


# # --- apply_variance_filter function ---
# def apply_variance_filter(df, n_top_features, data_type_name="Data"):
#     """Applies variance filtering to keep top N features."""
#     # Input df will be indexed by original_rna_sample_id after alignment
#     # Output df_filtered will also be indexed by original_rna_sample_id
#     print(f"\n--- Applying Variance Filter to {data_type_name} ---")
#     if df is None or df.empty: print("Input DataFrame empty."); return None
#     print(f"Input shape: {df.shape}. Targeting top {n_top_features} features.")
#     try:
#         df_numeric = df.select_dtypes(include=np.number)
#         if df_numeric.empty: print("Error: No numeric columns."); return None
#         variances = df_numeric.var(axis=0).dropna()
#         if variances.empty: print("Error: No variance calculated."); return None
#         n_to_keep = min(n_top_features, len(variances))
#         print(f"Selecting top {n_to_keep} features.")
#         top_features = variances.nlargest(n_to_keep).index.tolist()
#         df_filtered = df[top_features].copy()
#         print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
#         return df_filtered # Returns df indexed by original_rna_sample_id
#     except Exception as e: print(f"Error during variance filtering: {e}"); return None

# # --- REMOVED align_data_and_index function ---
# # --- REMOVED _format_id_component and _create_composite_id helpers ---
# # (As they are not needed with the original_rna_sample_id indexing strategy)