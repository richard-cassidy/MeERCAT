# mercat_analyzer/meercat/preprocess.py

import pandas as pd
import numpy as np
import re
# Import utils and config functions/variables
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
    missing_cols = [col for col in id_cols if col not in df.columns]
    if missing_cols:
        print(f"  ERROR: Missing required columns for composite ID: {missing_cols}")
        return None # Cannot create ID

    id_parts = []
    # Use provided experiment ID if available, otherwise try to get from column
    if exp_id is not None:
         id_parts.append(pd.Series(str(exp_id), index=df.index))
    elif 'experiment_id' in id_cols and 'experiment_id' in df.columns:
         id_parts.append(df['experiment_id'].astype(str).fillna('na'))
    else:
         print("  ERROR: Experiment ID missing for composite ID creation.")
         return None # Experiment ID is crucial

    # Format other components consistently
    for col in id_cols:
        if col != 'experiment_id': # Already handled exp_id
             id_parts.append(_format_id_component(df[col]))

    # Combine parts
    composite_id_series = pd.Series("_".join(map(str, tpl)) for tpl in zip(*id_parts))
    composite_id_series.index = df.index # Ensure index alignment
    print("  Composite ID Series created.")
    return composite_id_series

# --- Existing combine_dataframes function (Keep as before) ---
def combine_dataframes(df_dict, axis=0, join='outer'):
    # ... (function code remains the same) ...
    if not df_dict:
        print("Warning: No DataFrames provided for combination.")
        return None
    print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
    try:
        # Ensure keys are preserved if axis=1, handle potential integer keys
        if axis == 1:
             combined_df = pd.concat(df_dict, axis=axis, join=join)
        else: # axis=0
             combined_df = pd.concat(df_dict.values(), axis=axis, join=join, ignore_index=True)

        print(f"Combined shape: {combined_df.shape}")
        return combined_df
    except Exception as e:
        print(f"Error during concatenation: {e}")
        return None


# --- Revised clean_metabolite_data ---
def clean_metabolite_data(combined_df):
    """Cleans combined metabolite data: numeric conversion, NaN handling, filtering, indexing."""
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
    potential_metadata_cols = list(set(non_metabolite_cols + METADATA_COMPOSITE_ID_COLS)) # Include composite ID cols
    metadata_df = df_copy[[col for col in potential_metadata_cols if col in df_copy.columns]].copy()

    # 2. Convert features to numeric
    print("Converting features to numeric...")
    for col in metabolite_feature_cols:
        if col in df_copy.columns: # Check if column still exists
            if not pd.api.types.is_numeric_dtype(df_copy[col]):
                 df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce')
    # Keep only METADATA and NUMERIC metabolite feature columns
    numeric_feature_cols_present = df_copy[metabolite_feature_cols].select_dtypes(include=np.number).columns.tolist()
    cols_to_keep = list(metadata_df.columns) + numeric_feature_cols_present
    # Ensure no duplicates in cols_to_keep (metadata columns might overlap)
    cols_to_keep = sorted(list(set(cols_to_keep)))
    df_numeric_cleaned = df_copy[cols_to_keep].copy()
    print(f"Kept metadata and {len(numeric_feature_cols_present)} numeric feature columns.")


    # 3. Drop all-NaN columns (feature columns only)
    print("Checking feature columns for all-NaNs...")
    nan_feature_cols = [col for col in numeric_feature_cols_present if df_numeric_cleaned[col].isnull().all()]
    if nan_feature_cols:
        print(f"Identified {len(nan_feature_cols)} all-NaN feature columns to remove: {nan_feature_cols}")
        df_numeric_cleaned.drop(columns=nan_feature_cols, inplace=True)
        print(f"Dropped all-NaN feature columns. Shape now: {df_numeric_cleaned.shape}")
    else: print("No all-NaN feature columns found.")


    # 4. Filter out QC/Blank samples
    print("Filtering out QC/Blank samples...")
    df_filtered = df_numeric_cleaned # Start with current state
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
        return df_filtered, None # Return empty df, no metadata

    # Separate final metadata
    final_metadata_cols = [col for col in metadata_df.columns if col in df_filtered.columns]
    metadata_df_final = df_filtered[final_metadata_cols].copy()
    final_feature_cols = [col for col in df_filtered.columns if col not in final_metadata_cols]
    df_features_only = df_filtered[final_feature_cols].copy()


    # 5. Create and Set Composite Index (on both dataframes)
    print("Creating composite sample index...")
    index_set_successfully = False
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
            metadata_df_final.index = cleaned_final_index
            metadata_df_final.index.name = 'composite_sample_id'
            df_features_only.index = cleaned_final_index
            df_features_only.index.name = 'composite_sample_id'
            print(f"Successfully set index: '{df_features_only.index.name}'")
            index_set_successfully = True
    else:
         print("Failed to create composite ID series for metabolites.")


    print("Metabolite cleaning finished.")
    # Return the features DataFrame with the composite index, and the metadata DF
    if index_set_successfully:
         return df_features_only, metadata_df_final
    else:
         # Return features without index set, and metadata possibly with temp ID column
         return df_features_only, metadata_df_final


# --- Revised process_rna_dataframe ---
def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
    """Orients, cleans, labels, and indexes a single raw RNA DataFrame, separates metadata."""
    print(f"\nProcessing RNA DataFrame from '{filename}'...")
    if df_raw is None or df_raw.empty: return None, None

    metadata_df = None
    features_df = None

    try:
        # --- Basic Setup ---
        match_exp = re.match(r'(MC\d+)', filename, re.IGNORECASE) or re.search(r'(MC\d+)', filename, re.IGNORECASE)
        experiment_id = match_exp.group(1).upper() if match_exp else "UNKNOWN_EXPERIMENT"
        if experiment_id == "UNKNOWN_EXPERIMENT": print(f"  Warning: No MC ID found in filename '{filename}'.")

        # Assume first column is gene identifier if rows_are_genes, else index is gene
        if rows_are_genes:
            df_reset = df_raw.reset_index()
            gene_id_col = df_reset.columns[0]
            sample_cols = [col for col in df_reset.columns if col != gene_id_col]
            df_with_index = df_reset.set_index(gene_id_col)
            df_numeric_part = df_with_index[sample_cols]
            features_df = df_numeric_part.T.copy() # Samples x Genes
        else:
            features_df = df_raw.copy() # Assumes Samples x Genes already
        print(f"  Oriented shape (Samples x Features): {features_df.shape}")

        # Clean Index (Sample IDs) and Columns (Gene IDs)
        features_df.index = clean_index(features_df.index)
        features_df.columns = clean_index(features_df.columns)

        # Clean Gene IDs & Handle Duplicates
        features_df.columns = clean_feature_names(features_df.columns)
        features_df = handle_duplicate_features(features_df) # Returns only numeric columns

        # --- Extract Metadata (Simplified - Requires External Metadata for Robustness) ---
        # Create a basic metadata df with just the experiment ID and original index
        metadata_df = pd.DataFrame(index=features_df.index)
        metadata_df['experiment_id'] = experiment_id
        metadata_df['original_rna_sample_id'] = features_df.index # Store for reference

        # **** Attempt to parse index IF NEEDED (Less robust) ****
        # print("  Attempting metadata extraction from index (heuristic)...")
        # Add your regex parsing here if you MUST rely on it, creating columns
        # in metadata_df like 'replicate', 'arsenic_concentration', 'days'.
        # Ensure these columns match METADATA_COMPOSITE_ID_COLS expected names ('ias_conc', 'day', 'rep')
        # If parsing fails, these columns might be missing or NaN.
        # Example (replace with your actual parsing):
        # pattern = ... ; extracted_data = metadata_df['original_rna_sample_id'].str.extract(...)
        # metadata_df['rep'] = pd.to_numeric(extracted_data['rep_col_name'], errors='coerce')
        # metadata_df['ias_conc'] = pd.to_numeric(extracted_data['conc_col_name'], errors='coerce')
        # metadata_df['day'] = pd.to_numeric(extracted_data['day_col_name'], errors='coerce')
        print("  WARNING: RNA Metadata extraction from index is not implemented robustly.")
        print("           Composite index relies only on info available (e.g., experiment_id).")
        print("           RECOMMEND providing full metadata externally for RNA samples.")
        # Create placeholder columns if parsing not done/failed, to match expected ID cols
        for col in METADATA_COMPOSITE_ID_COLS:
             if col not in metadata_df.columns and col != 'experiment_id':
                  metadata_df[col] = np.nan


        # --- Create and Set Composite Index ---
        print("  Creating composite sample ID...")
        index_set_successfully = False
        composite_id_series = _create_composite_id(metadata_df, METADATA_COMPOSITE_ID_COLS, exp_id=experiment_id) # Pass exp_id

        if composite_id_series is not None:
            if composite_id_series.duplicated().any():
                print(f"  WARNING: Duplicate composite IDs generated for RNA! Index not set.")
                metadata_df['composite_sample_id_TEMP'] = composite_id_series
                features_df['composite_sample_id_TEMP'] = composite_id_series # Add here too if index fails
            else:
                print("  Composite ID unique. Setting as index...")
                cleaned_final_index = clean_index(composite_id_series)
                metadata_df.index = cleaned_final_index
                metadata_df.index.name = 'composite_sample_id'
                features_df.index = cleaned_final_index
                features_df.index.name = 'composite_sample_id'
                print(f"  Set composite ID as index. Example: {features_df.index[0]}")
                index_set_successfully = True
        else:
            print("  Failed to create composite ID series for RNA.")

        # Convert features to numeric (do this *after* potential handling of duplicates)
        features_df = features_df.apply(pd.to_numeric, errors='coerce')

        print(f"Finished processing RNA file {filename}.")
        if index_set_successfully:
             return features_df, metadata_df
        else: # Return un-indexed if ID creation failed
             print("  Returning RNA data without composite index.")
             return features_df, metadata_df


    except Exception as e:
         print(f"  An unexpected error occurred during RNA processing for '{filename}': {e}")
         return None, None

# --- Revised normalize_rna_data ---
def normalize_rna_data(combined_processed_rna):
    """Filters low-count genes and performs CPM + log2(x+1) normalization. Assumes input has features only."""
    print("\n--- Filtering and Normalizing Combined RNA Data ---")
    if combined_processed_rna is None or combined_processed_rna.empty:
        print("Input DataFrame is None or empty. Skipping.")
        return None

    # Assume input is features only (numeric)
    rna_counts = combined_processed_rna.select_dtypes(include=np.number).fillna(0).astype(np.float32)
    if rna_counts.empty:
        print("Error: No numeric gene columns identified in combined RNA data.")
        return None
    print(f"Input count data shape (Samples x Genes): {rna_counts.shape}")

    # 1. Filtering
    min_samples = max(2, int(len(rna_counts) * RNA_FILTER_MIN_SAMPLES_FRAC))
    print(f"Filtering: Keeping genes with >= {RNA_FILTER_MIN_COUNTS} counts in at least {min_samples} samples.")
    counts_array = rna_counts.values
    genes_to_keep_mask = (counts_array >= RNA_FILTER_MIN_COUNTS).sum(axis=0) >= min_samples
    rna_filtered = rna_counts.loc[:, genes_to_keep_mask]
    print(f"Removed {rna_counts.shape[1] - rna_filtered.shape[1]} genes. Shape after filtering: {rna_filtered.shape}")

    if rna_filtered.empty:
        print("Error: Filtering removed all genes.")
        return None

    # 2. Normalization (CPM + Log2)
    print("Normalizing filtered counts (Vectorized CPM + Log2)...")
    counts_filtered_array = rna_filtered.values
    library_sizes = counts_filtered_array.sum(axis=1)
    library_sizes_col = library_sizes[:, np.newaxis] + 1e-9 # Add epsilon for stability
    cpm_array = (counts_filtered_array / library_sizes_col) * 1e6
    log2_cpm_array = np.log2(cpm_array + 1)

    rna_normalized_log2 = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
    print("Normalization complete.")
    print("Final normalized RNA data shape:", rna_normalized_log2.shape)
    return rna_normalized_log2


# --- Existing align_samples (Keep as before, relies on index) ---
def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
    # ... (function code remains the same) ...
    print(f"\n--- Aligning {df1_name} and {df2_name} to Common Samples ---")
    if df1 is None or df2 is None or df1.empty or df2.empty:
        print("One or both input DataFrames are missing or empty. Cannot align.")
        return None, None
    if not isinstance(df1.index, pd.Index) or not isinstance(df2.index, pd.Index):
         print("Error: One or both DataFrames lack a valid index.")
         return None, None

    print(f"Input {df1_name} shape: {df1.shape}")
    print(f"Input {df2_name} shape: {df2.shape}")

    # Clean indices before intersection
    index1_clean = clean_index(df1.index)
    index2_clean = clean_index(df2.index)
    common_samples = index1_clean.intersection(index2_clean)
    print(f"Found {len(common_samples)} common samples based on cleaned index.")

    if len(common_samples) == 0:
        print("CRITICAL WARNING: No common samples found!")
        return None, None

    try:
        # Filter using original index corresponding to cleaned common samples
        # Need to ensure the indices used here ARE the cleaned ones for lookup if df1/df2 haven't been reindexed yet
        df1_matched = df1[index1_clean.isin(common_samples)].copy()
        df2_matched = df2[index2_clean.isin(common_samples)].copy()

        # Re-clean indices after filtering and align
        df1_matched.index = clean_index(df1_matched.index)
        df2_matched.index = clean_index(df2_matched.index)

        # Align df2 to df1's index order
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
    # ... (function code remains the same) ...
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
             # Return original df? Or empty? Let's return None for consistency
             return None

        variances = df_numeric.var(axis=0).dropna() # Drop NaNs if variance is NaN
        if variances.empty:
             print("Error: Could not calculate variance for any feature.")
             return None # Return None

        n_to_keep = min(n_top_features, len(variances))
        print(f"Selecting top {n_to_keep} features (out of {len(variances)} with calculated variance).")
        top_features = variances.nlargest(n_to_keep).index.tolist()

        # Filter the original DataFrame using the selected feature names
        # This preserves any non-numeric columns if they exist
        df_filtered = df[top_features].copy()
        print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
        return df_filtered

    except Exception as e:
        print(f"Error during variance filtering for {data_type_name}: {e}")
        return None