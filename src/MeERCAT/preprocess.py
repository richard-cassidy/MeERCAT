# mercat_analyzer/meercat/preprocess.py

import pandas as pd
import numpy as np
import re
# Import utils and config functions/variables
# Assuming utils.py is in the same directory ('.') and config.py as well
from .utils import identify_metabolite_columns, clean_feature_names, handle_duplicate_features, clean_index
from .config import (
    METADATA_SAMPLE_COL, SAMPLE_FILTER_EXCLUDE_PATTERNS,
    METADATA_COMPOSITE_ID_COLS, RNA_FILTER_MIN_COUNTS,
    RNA_FILTER_MIN_SAMPLES_FRAC, VAR_FILTER_TOP_N_GENES,
    VAR_FILTER_TOP_N_METABOLITES, METABOLITE_ID_PREFIXES,
    METABOLITE_ID_SUFFIX_REGEX, RNA_METADATA_COLS_EXPECTED
)

# --- Helper for Consistent ID String Formatting ---
def _format_id_component(series):
    """Converts a numeric series to string for IDs, handling NaN and '.0'."""
    num_series = pd.to_numeric(series, errors='coerce')
    # Apply formatting: int if whole number, else string, 'na' if NaN
    str_series = num_series.apply(
        lambda x: str(int(x)) if pd.notna(x) and x == int(x) else str(x) if pd.notna(x) else 'na'
    )
    return str_series

# --- Helper to Create Composite ID ---
def _create_composite_id(df, id_cols, exp_id=None):
    """Creates a composite ID Series from specified columns in a DataFrame."""
    # id_cols expected: ['experiment_id', 'ias_conc', 'day', 'rep'] from config
    print(f"  Attempting to create composite ID using columns: {id_cols}")
    id_parts = []

    # --- Handle Experiment ID ---
    # Use provided experiment ID if available (e.g., from filename)
    if exp_id is not None:
         id_parts.append(pd.Series(str(exp_id), index=df.index))
         print(f"    Using provided experiment_id: {exp_id}")
         # Remove 'experiment_id' from id_cols if it exists, as we used the override
         id_cols = [col for col in id_cols if col != 'experiment_id']
    # Otherwise, check if 'experiment_id' column exists in df
    elif 'experiment_id' in id_cols and 'experiment_id' in df.columns:
         id_parts.append(df['experiment_id'].astype(str).fillna('na'))
         print("    Using 'experiment_id' column from DataFrame.")
         # Remove 'experiment_id' from id_cols to avoid double processing
         id_cols = [col for col in id_cols if col != 'experiment_id']
    else:
         print("  ERROR: Experiment ID column ('experiment_id') not found in DataFrame and not provided separately.")
         return None # Experiment ID is crucial

    # --- Check for remaining required columns ---
    missing_cols = [col for col in id_cols if col not in df.columns]
    if missing_cols:
        print(f"  ERROR: Missing required columns for composite ID: {missing_cols}")
        return None # Cannot create ID

    # --- Format other components consistently ---
    print(f"    Processing components: {id_cols}")
    for col in id_cols:
        id_parts.append(_format_id_component(df[col]))

    # --- Combine parts ---
    try:
        composite_id_series = pd.Series("_".join(map(str, tpl)) for tpl in zip(*id_parts))
        composite_id_series.index = df.index # Ensure index alignment
        print("  Composite ID Series created successfully.")
        return composite_id_series
    except Exception as e_zip:
         print(f"  ERROR during ID component combination: {e_zip}")
         return None


# --- Existing combine_dataframes function (Keep as before) ---
def combine_dataframes(df_dict, axis=0, join='outer'):
    """Combines a dictionary of DataFrames."""
    if not df_dict:
        print("Warning: No DataFrames provided for combination.")
        return None
    print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
    try:
        # Ensure keys are preserved if axis=1, handle potential integer keys
        if axis == 1:
             combined_df = pd.concat(df_dict, axis=axis, join=join)
        else: # axis=0
             # Check if values are DataFrames before concatenating
             valid_dfs = [df for df in df_dict.values() if isinstance(df, pd.DataFrame)]
             if not valid_dfs:
                  print("Warning: No valid DataFrames found in dictionary values.")
                  return None
             combined_df = pd.concat(valid_dfs, axis=axis, join=join, ignore_index=(axis==0))

        print(f"Combined shape: {combined_df.shape}")
        return combined_df
    except Exception as e:
        print(f"Error during concatenation: {e}")
        return None


# --- Revised clean_metabolite_data ---
def clean_metabolite_data(combined_df):
    """Cleans combined metabolite data: separates metadata, cleans features, filters, creates index.
       Returns:
           df_features_only (DataFrame): Numeric features indexed by composite_sample_id.
           metadata_df_final (DataFrame): Metadata indexed by composite_sample_id.
    """
    print("\n--- Cleaning Combined Metabolite Data ---")
    if combined_df is None or combined_df.empty:
        print("Input DataFrame is None or empty. Skipping cleaning.")
        return None, None # Return None for df and metadata

    df_copy = combined_df.copy()
    metadata_df = None # Initialize metadata df

    # 1. Identify columns
    all_cols = df_copy.columns.tolist()
    metabolite_feature_cols, non_metabolite_cols = identify_metabolite_columns(
        all_cols, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX
    )
    print(f"Identified {len(metabolite_feature_cols)} potential metabolite features.")
    print(f"Found {len(non_metabolite_cols)} assumed non-feature columns (metadata/IDs): {sorted([str(c) for c in non_metabolite_cols])}")

    # Extract potential metadata BEFORE cleaning features
    # Include columns needed for composite ID creation explicitly
    potential_metadata_cols = list(set(non_metabolite_cols + METADATA_COMPOSITE_ID_COLS + [METADATA_SAMPLE_COL]))
    metadata_cols_present = [col for col in potential_metadata_cols if col in df_copy.columns]
    if metadata_cols_present:
        metadata_df = df_copy[metadata_cols_present].copy()
        print(f"Extracted potential metadata columns: {metadata_cols_present}")
    else:
        print("Warning: No potential metadata columns found based on config.")
        metadata_df = pd.DataFrame(index=df_copy.index) # Create empty metadata df with index

    # 2. Convert features to numeric
    print("Converting features to numeric...")
    feature_cols_present = [col for col in metabolite_feature_cols if col in df_copy.columns]
    df_features = df_copy[feature_cols_present].copy() # Work on features separately
    converted_count = 0
    for col in feature_cols_present:
        if not pd.api.types.is_numeric_dtype(df_features[col]):
            original_dtype = df_features[col].dtype
            df_features[col] = pd.to_numeric(df_features[col], errors='coerce')
            if original_dtype != df_features[col].dtype: converted_count += 1
    print(f"Numeric conversion attempted ({converted_count} columns potentially changed).")

    # Select only numeric features after conversion
    numeric_feature_cols_present = df_features.select_dtypes(include=np.number).columns.tolist()
    df_numeric_features = df_features[numeric_feature_cols_present].copy()
    print(f"Kept {len(numeric_feature_cols_present)} numeric feature columns.")


    # 3. Drop all-NaN columns (feature columns only)
    print("Checking feature columns for all-NaNs...")
    nan_feature_cols = [col for col in numeric_feature_cols_present if df_numeric_features[col].isnull().all()]
    if nan_feature_cols:
        print(f"Identified {len(nan_feature_cols)} all-NaN feature columns to remove: {nan_feature_cols}")
        df_numeric_features.drop(columns=nan_feature_cols, inplace=True)
        print(f"Dropped all-NaN feature columns. Features shape now: {df_numeric_features.shape}")
    else: print("No all-NaN feature columns found.")

    # Re-attach metadata to the cleaned numeric features for sample filtering
    df_to_filter = pd.merge(metadata_df, df_numeric_features, left_index=True, right_index=True, how='inner')
    print(f"Re-merged metadata and cleaned numeric features. Shape: {df_to_filter.shape}")

    # 4. Filter out QC/Blank samples
    print("Filtering out QC/Blank samples...")
    df_filtered = df_to_filter # Start with current state
    if METADATA_SAMPLE_COL in df_filtered.columns:
        exclude_regex = '|'.join(SAMPLE_FILTER_EXCLUDE_PATTERNS)
        sample_col_str = df_filtered[METADATA_SAMPLE_COL].astype(str)
        keep_mask = ~sample_col_str.str.lower().str.contains(exclude_regex, na=False, regex=True)
        rows_before = len(df_filtered)
        df_filtered = df_filtered[keep_mask].copy()
        print(f"Filtered {rows_before - len(df_filtered)} QC/Blank/Other rows. Shape now: {df_filtered.shape}")
    else: print(f"Warning: '{METADATA_SAMPLE_COL}' column not found for filtering.")

    if df_filtered.empty:
        print("Warning: DataFrame empty after filtering samples.")
        return df_filtered, pd.DataFrame() # Return empty dfs

    # 5. Separate final metadata and features AGAIN after sample filtering
    final_metadata_cols = [col for col in metadata_cols_present if col in df_filtered.columns]
    metadata_df_final = df_filtered[final_metadata_cols].copy()
    final_feature_cols = [col for col in df_filtered.columns if col not in final_metadata_cols]
    df_features_only = df_filtered[final_feature_cols].copy()
    print(f"Separated final Metadata ({metadata_df_final.shape}) and Features ({df_features_only.shape})")

    # 6. Create and Set Composite Index (on both dataframes)
    print("Creating composite sample index...")
    index_set_successfully = False
    # Use metadata_df_final which should contain the necessary columns
    composite_id_series = _create_composite_id(metadata_df_final, METADATA_COMPOSITE_ID_COLS)

    if composite_id_series is not None:
        if composite_id_series.duplicated().any():
            print("CRITICAL WARNING: Duplicate composite IDs generated for metabolites! Index cannot be set reliably.")
            # Optionally store ID as a column for inspection
            metadata_df_final['composite_sample_id_TEMP'] = composite_id_series
            df_features_only['composite_sample_id_TEMP'] = composite_id_series
        else:
            print("Composite ID is unique. Setting as index...")
            cleaned_final_index = clean_index(composite_id_series)
            try:
                 metadata_df_final.index = cleaned_final_index
                 metadata_df_final.index.name = 'composite_sample_id'
                 df_features_only.index = cleaned_final_index
                 df_features_only.index.name = 'composite_sample_id'
                 print(f"Successfully set index: '{df_features_only.index.name}'")
                 index_set_successfully = True
            except Exception as e_index:
                 print(f"ERROR setting index after creation: {e_index}")
    else:
         print("Failed to create composite ID series for metabolites.")


    print("Metabolite cleaning finished.")
    # Return the features DataFrame and the metadata DF, with index potentially set
    if index_set_successfully:
         return df_features_only, metadata_df_final
    else:
         # Return features/metadata potentially without index set, may have TEMP ID column
         print("Warning: Returning data without composite index set due to errors/duplicates.")
         return df_features_only, metadata_df_final


# --- Revised process_rna_dataframe ---
def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
    """Orients, cleans features, and returns RNA features indexed by original sample ID.
       Assumes metadata will be merged externally.
       Returns:
           features_df (DataFrame): Numeric features indexed by original_rna_sample_id.
    """
    print(f"\nProcessing RNA DataFrame from '{filename}'...")
    if df_raw is None or df_raw.empty: return None

    features_df = None

    try:
        # --- Basic Setup & Orientation ---
        if rows_are_genes:
            df_reset = df_raw.reset_index()
            gene_id_col = df_reset.columns[0]
            sample_cols = [col for col in df_reset.columns if col != gene_id_col]
            if not sample_cols: raise ValueError("No sample columns found.")
            df_with_index = df_reset.set_index(gene_id_col)
            # Convert to numeric EARLIER, handle errors
            df_numeric_part = df_with_index[sample_cols].apply(pd.to_numeric, errors='coerce')
            features_df = df_numeric_part.T.copy() # Samples x Genes
        else:
            # Assumes Samples x Genes already, convert to numeric
            features_df = df_raw.apply(pd.to_numeric, errors='coerce').copy()

        print(f"  Oriented shape (Samples x Features): {features_df.shape}")

        # --- Clean Index (Original Sample IDs from file) ---
        features_df.index = clean_index(features_df.index)
        features_df.index.name = 'original_rna_sample_id' # Give it a meaningful name

        # --- Clean Columns (Gene IDs) & Handle Duplicates ---
        features_df.columns = clean_index(features_df.columns) # Clean gene IDs too
        features_df.columns = clean_feature_names(features_df.columns)
        features_df = handle_duplicate_features(features_df) # Returns only numeric columns

        # --- Fill NaNs introduced by numeric conversion or duplicates ---
        # (handle_duplicate_features already fills with 0 after sum)
        # If NaNs could remain from initial pd.to_numeric, fill them here
        nan_count = features_df.isna().sum().sum()
        if nan_count > 0:
             print(f"  Filling {nan_count} NaNs in features with 0.")
             features_df.fillna(0, inplace=True)

        print(f"Finished basic processing for RNA file {filename}. Index is original sample ID.")
        # Composite index creation and metadata merge happens LATER in the main script
        return features_df # Return only the features dataframe

    except Exception as e:
         print(f"  An unexpected error occurred during basic RNA processing for '{filename}': {e}")
         return None


# --- Revised normalize_rna_data ---
def normalize_rna_data(rna_features_indexed):
    """Filters low-count genes and performs CPM + log2(x+1) normalization.
       Assumes input is features only, indexed by composite_sample_id.
    """
    print("\n--- Filtering and Normalizing Combined RNA Data ---")
    if rna_features_indexed is None or rna_features_indexed.empty:
        print("Input DataFrame is None or empty. Skipping normalization.")
        return None

    # Assume input is features only (numeric checked before this step ideally)
    rna_counts = rna_features_indexed.select_dtypes(include=np.number).fillna(0).astype(np.float32)
    if rna_counts.empty:
        print("Error: No numeric gene columns identified for normalization.")
        return None
    print(f"Input count data shape (Samples x Genes): {rna_counts.shape}")

    # 1. Filtering
    min_samples = max(2, int(len(rna_counts) * RNA_FILTER_MIN_SAMPLES_FRAC))
    print(f"Filtering: Keeping genes with >= {RNA_FILTER_MIN_COUNTS} counts in at least {min_samples} samples.")
    try:
        counts_array = rna_counts.values
        genes_to_keep_mask = (counts_array >= RNA_FILTER_MIN_COUNTS).sum(axis=0) >= min_samples
        rna_filtered = rna_counts.loc[:, genes_to_keep_mask]
        print(f"Removed {rna_counts.shape[1] - rna_filtered.shape[1]} genes. Shape after filtering: {rna_filtered.shape}")
    except Exception as e_filt_mask:
         print(f"Error during gene filtering mask creation: {e_filt_mask}")
         return None

    if rna_filtered.empty:
        print("Error: Filtering removed all genes.")
        return None

    # 2. Normalization (CPM + Log2)
    print("Normalizing filtered counts (Vectorized CPM + Log2)...")
    try:
        counts_filtered_array = rna_filtered.values
        library_sizes = counts_filtered_array.sum(axis=1)
        # Add epsilon for stability BEFORE division
        library_sizes_col = library_sizes[:, np.newaxis] + 1e-9
        # Check for zero library sizes after adding epsilon (shouldn't happen)
        if np.any(library_sizes_col <= 1e-9):
             print("Warning: Some library sizes are zero or near-zero after adding epsilon!")
             # Handle this? e.g., replace those rows/cols with zeros?
             # For now, proceed, log2(1) will be 0 for these.

        cpm_array = (counts_filtered_array / library_sizes_col) * 1e6
        # Ensure non-negativity before log (should be fine with +1)
        log2_cpm_array = np.log2(np.maximum(cpm_array, 0) + 1) # Use np.maximum for safety

        rna_normalized_log2 = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
        print("Normalization complete.")
        print("Final normalized RNA data shape:", rna_normalized_log2.shape)
        return rna_normalized_log2
    except Exception as e_norm:
        print(f"Error during normalization: {e_norm}")
        return None


# --- Revised align_samples ---
def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
    """Aligns two DataFrames based on common indices. Assumes indices are composite_sample_id."""
    print(f"\n--- Aligning {df1_name} and {df2_name} to Common Samples ---")
    if df1 is None or df2 is None or df1.empty or df2.empty:
        print("One or both input DataFrames are missing or empty. Cannot align.")
        return None, None
    # Check if indices are already set (should be composite IDs here)
    if df1.index.name != 'composite_sample_id' or df2.index.name != 'composite_sample_id':
         print("Error: Input DataFrames for alignment must be indexed by 'composite_sample_id'.")
         print(f"       Index names found: '{df1.index.name}', '{df2.index.name}'")
         return None, None

    print(f"Input {df1_name} shape: {df1.shape}")
    print(f"Input {df2_name} shape: {df2.shape}")

    # Indices should already be clean composite IDs, but run clean_index just in case
    index1_clean = clean_index(df1.index)
    index2_clean = clean_index(df2.index)
    common_samples = index1_clean.intersection(index2_clean)
    print(f"Found {len(common_samples)} common samples based on composite index.")

    if len(common_samples) == 0:
        print("CRITICAL WARNING: No common samples found based on composite index!")
        return None, None

    try:
        # Filter using the common index directly
        df1_matched = df1.loc[common_samples].copy()
        df2_matched = df2.loc[common_samples].copy()

        # Reindex df2 to match df1's exact order
        df2_aligned = df2_matched.reindex(df1_matched.index)

        # Final verification
        if not df1_matched.index.equals(df2_aligned.index):
             raise ValueError("Final indices do not match after alignment!")

        print(f"Final aligned {df1_name} shape: {df1_matched.shape}")
        print(f"Final aligned {df2_name} shape: {df2_aligned.shape}")
        print("Alignment successful.")
        return df1_matched, df2_aligned

    except Exception as e:
        print(f"Error during sample alignment: {e}")
        return None, None


# --- Existing apply_variance_filter (Keep as before) ---
def apply_variance_filter(df, n_top_features, data_type_name="Data"):
    """Applies variance filtering to keep top N features."""
    print(f"\n--- Applying Variance Filter to {data_type_name} ---")
    if df is None or df.empty:
        print("Input DataFrame is None or empty. Skipping.")
        return None

    print(f"Input shape: {df.shape}")
    print(f"Targeting top {n_top_features} features.")

    try:
        # Select numeric features only for variance calculation
        df_numeric = df.select_dtypes(include=np.number)
        if df_numeric.empty:
             print("Error: No numeric columns found for variance filtering.")
             return None # Return None if no numeric data

        variances = df_numeric.var(axis=0).dropna() # Drop NaNs if variance is NaN
        if variances.empty:
             print("Error: Could not calculate variance for any feature.")
             return None # Return None

        n_to_keep = min(n_top_features, len(variances))
        print(f"Selecting top {n_to_keep} features (out of {len(variances)} with calculated variance).")
        top_features = variances.nlargest(n_to_keep).index.tolist()

        # Filter the DataFrame using the selected feature names
        df_filtered = df[top_features].copy() # Use copy
        print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
        return df_filtered

    except Exception as e:
        print(f"Error during variance filtering for {data_type_name}: {e}")
        return None