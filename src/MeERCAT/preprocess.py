# mercat_analyzer/meercat/preprocess.py

import pandas as pd
import numpy as np
import re
# Import utils and config functions/variables
from .utils import identify_metabolite_columns, clean_feature_names, handle_duplicate_features
from .config import (
    METADATA_SAMPLE_COL, SAMPLE_FILTER_EXCLUDE_PATTERNS,
    # METADATA_COMPOSITE_ID_COLS, # No longer needed for index creation here
    RNA_FILTER_MIN_COUNTS,
    RNA_FILTER_MIN_SAMPLES_FRAC, VAR_FILTER_TOP_N_GENES,
    VAR_FILTER_TOP_N_METABOLITES, METABOLITE_ID_PREFIXES,
    METABOLITE_ID_SUFFIX_REGEX #, RNA_METADATA_COLS_EXPECTED # Not used directly here now
)
# Import clean_index from load_data
try:
    from .load_data import clean_index
except ImportError:
    print("ERROR in preprocess.py: Could not import clean_index from .load_data.")
    def clean_index(index): print("ERROR: clean_index function missing!"); return index


# --- combine_dataframes function ---
def combine_dataframes(df_dict, axis=0, join='outer'):
    """Combines a dictionary of DataFrames."""
    if not isinstance(df_dict, dict) or not df_dict: print("Warning: Input not valid dict."); return None
    print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
    try:
        valid_dfs = [df for df in df_dict.values() if isinstance(df, pd.DataFrame) and not df.empty]
        if not valid_dfs: print("Warning: No valid DataFrames found."); return None
        # Preserve index when combining rows (axis=0)
        combined_df = pd.concat(valid_dfs, axis=axis, join=join, ignore_index=False) # ignore_index=False
        print(f"Combined shape: {combined_df.shape}")
        # Use INFO for duplicate index warning as it's expected when combining raw files
        if axis == 0 and combined_df.index.has_duplicates: print(f"INFO: Input indices had duplicates (expected for axis=0 combine).")
        return combined_df
    except Exception as e: print(f"Error during concatenation: {e}"); return None


# --- clean_metabolite_data function ---
def clean_metabolite_data(combined_df):
    """Basic cleaning: averages duplicate samples (based on 'original_rna_sample_id'),
       separates metadata, cleans features, filters samples,
       and sets index using 'original_rna_sample_id'.
       Returns:
           df_features_indexed (DataFrame): Numeric features indexed by original_rna_sample_id.
           metadata_df_indexed (DataFrame): Metadata indexed by original_rna_sample_id.
    """
    print("\n--- Cleaning Combined Metabolite Data ---")
    if combined_df is None or combined_df.empty: return None, None
    df_cleaned = combined_df.copy(); metadata_df_final = None; df_features_only = None; index_set_successfully = False

    INDEX_COL_NAME = 'original_rna_sample_id' # This column MUST exist in combined_df

    try:
        if INDEX_COL_NAME not in df_cleaned.columns:
             raise ValueError(f"Crucial column '{INDEX_COL_NAME}' not found in combined metabolite data.")

        # --- Handle Duplicates by Averaging Features ---
        print(f"Checking for duplicates in '{INDEX_COL_NAME}' column...")
        df_cleaned[INDEX_COL_NAME] = clean_index(df_cleaned[INDEX_COL_NAME]) # Clean the ID column first
        duplicates_present = df_cleaned[INDEX_COL_NAME].duplicated().any()
        index_is_unique = False # Assume not unique initially if duplicates found

        if duplicates_present:
            print(f"  Duplicates found in '{INDEX_COL_NAME}'. Averaging numeric features...")
            # Identify numeric feature columns BEFORE grouping
            feature_cols_initial, non_feature_cols_initial = identify_metabolite_columns(df_cleaned.columns, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX)
            numeric_feature_cols = df_cleaned[feature_cols_initial].select_dtypes(include=np.number).columns.tolist()
            # Grouping columns: The unique ID + any other metadata to keep 'first'
            metadata_cols_to_keep = [col for col in non_feature_cols_initial if col != INDEX_COL_NAME and col in df_cleaned.columns]
            # Make sure metadata columns are strings before 'first' aggregation if needed
            for col in metadata_cols_to_keep:
                 if not pd.api.types.is_numeric_dtype(df_cleaned[col]):
                      df_cleaned[col] = df_cleaned[col].astype(str)

            agg_dict = {**{col: 'first' for col in metadata_cols_to_keep},
                        **{col: 'mean' for col in numeric_feature_cols}}

            # Perform grouping
            df_grouped = df_cleaned.groupby(INDEX_COL_NAME).agg(agg_dict)
            # Check if grouping resulted in an empty dataframe (unlikely but possible)
            if df_grouped.empty: raise ValueError("Grouping by ID resulted in an empty DataFrame.")
            # Check for duplicates in the NEW index (should not happen after groupby)
            if df_grouped.index.duplicated().any(): raise ValueError("Duplicates STILL EXIST after averaging! Check grouping columns/logic.")

            df_cleaned = df_grouped.reset_index() # Bring ID back as column for rest of processing
            print(f"  Shape after averaging duplicates: {df_cleaned.shape}")
            index_is_unique = True # Should be unique now
        else:
            print(f"  Values in '{INDEX_COL_NAME}' column are unique.")
            index_is_unique = True # Already unique

        # --- Continue with cleaning on the (potentially averaged) df_cleaned ---
        all_cols = df_cleaned.columns.tolist()
        metabolite_feature_cols, non_metabolite_cols = identify_metabolite_columns(all_cols, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX)
        metadata_cols_present = sorted(list(set([col for col in non_metabolite_cols if col != INDEX_COL_NAME] + [METADATA_SAMPLE_COL])))
        metadata_cols_present = [col for col in metadata_cols_present if col in df_cleaned.columns] # Ensure they exist

        print("Converting features to numeric (post-averaging check)...")
        feature_cols_present = [col for col in metabolite_feature_cols if col in df_cleaned.columns]
        if not feature_cols_present: print("Warning: No metabolite features remaining."); df_numeric_features = pd.DataFrame(index=df_cleaned.index)
        else:
            df_features = df_cleaned[feature_cols_present].copy()
            for col in feature_cols_present: # Convert again if averaging created non-numeric
                if not pd.api.types.is_numeric_dtype(df_features[col]): df_features[col] = pd.to_numeric(df_features[col], errors='coerce')
            numeric_feature_cols = df_features.select_dtypes(include=np.number).columns.tolist()
            df_numeric_features = df_features[numeric_feature_cols].copy()
            print(f"Kept {len(numeric_feature_cols)} numeric feature columns.")
            nan_feature_cols = [col for col in numeric_feature_cols if df_numeric_features[col].isnull().all()]
            if nan_feature_cols: df_numeric_features.drop(columns=nan_feature_cols, inplace=True); print(f"Dropped {len(nan_feature_cols)} all-NaN features.")

        # Prepare for sample filtering (re-attach relevant metadata)
        # Use the ID column to align indices before merging
        cols_for_meta_df = [INDEX_COL_NAME] + metadata_cols_present
        metadata_df = df_cleaned[cols_for_meta_df].set_index(INDEX_COL_NAME) # Index meta by ID
        df_numeric_features[INDEX_COL_NAME] = df_cleaned[INDEX_COL_NAME] # Ensure ID col exists
        df_numeric_features = df_numeric_features.set_index(INDEX_COL_NAME) # Index features by ID
        # Use inner join to ensure alignment and drop samples missing in either part after averaging/cleaning
        df_to_filter = pd.merge(metadata_df, df_numeric_features, left_index=True, right_index=True, how='inner')
        print(f"Re-merged metadata and numeric features. Shape: {df_to_filter.shape}")


        print("Filtering out QC/Blank samples...")
        df_filtered = df_to_filter
        if METADATA_SAMPLE_COL in df_filtered.columns:
            exclude_regex = '|'.join(SAMPLE_FILTER_EXCLUDE_PATTERNS)
            keep_mask = ~df_filtered[METADATA_SAMPLE_COL].astype(str).str.lower().str.contains(exclude_regex, na=False, regex=True)
            rows_before = len(df_filtered)
            df_filtered = df_filtered[keep_mask].copy()
            print(f"Filtered {rows_before - len(df_filtered)} QC/Blank rows. Shape now: {df_filtered.shape}")
        if df_filtered.empty: print("Warning: DataFrame empty after filtering."); return pd.DataFrame(), pd.DataFrame()

        # Final separation (index is now original_rna_sample_id)
        final_metadata_cols = [col for col in metadata_cols_present if col in df_filtered.columns]
        metadata_df_final = df_filtered[final_metadata_cols].copy()
        final_feature_cols = [col for col in df_filtered.columns if col not in final_metadata_cols]
        df_features_only = df_filtered[final_feature_cols].copy()
        print(f"Separated final Metadata ({metadata_df_final.shape}) and Features ({df_features_only.shape})")

        # Confirm index name
        df_features_only.index.name = INDEX_COL_NAME
        metadata_df_final.index.name = INDEX_COL_NAME
        index_set_successfully = True
        print(f"Index '{INDEX_COL_NAME}' confirmed on output DataFrames.")

    except Exception as e_clean: print(f"ERROR during metabolite cleaning: {e_clean}"); return None, None

    print("Metabolite cleaning finished.")
    # Return features and metadata indexed by original_rna_sample_id
    return df_features_only, metadata_df_final


# --- process_rna_dataframe function ---
def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
    """Orients, cleans features, returns RNA features indexed by original sample ID."""
    print(f"\nProcessing RNA DataFrame from '{filename}'...")
    if df_raw is None or df_raw.empty: return None
    features_df = None
    try:
        if rows_are_genes:
            df_reset = df_raw.reset_index(); gene_id_col = df_reset.columns[0]
            sample_pattern = r'^\d+_\d{2,}_\d+[a-zA-Z]+$'
            sample_cols = [col for col in df_reset.columns if col != gene_id_col and isinstance(col, str) and re.match(sample_pattern, col)]
            annotation_cols = [col for col in df_reset.columns if col != gene_id_col and col not in sample_cols]
            if not sample_cols: sample_cols = [col for col in df_reset.columns if col != gene_id_col]; print(f"  Warning: Using fallback for RNA sample columns in {filename}.")
            # if annotation_cols: print(f"  Identified {len(annotation_cols)} annotation columns: {annotation_cols}") # Less verbose
            df_numeric_part = df_reset.set_index(gene_id_col)[sample_cols].apply(pd.to_numeric, errors='coerce')
            features_df = df_numeric_part.T.copy()
        else: features_df = df_raw.select_dtypes(include=np.number).apply(pd.to_numeric, errors='coerce').copy()
        print(f"  Oriented features shape (Samples x Genes): {features_df.shape}")
        features_df.index = clean_index(features_df.index) # Uses imported clean_index
        features_df.index.name = 'original_rna_sample_id' # Set index name
        features_df.columns = clean_index(features_df.columns); features_df.columns = clean_feature_names(features_df.columns); features_df = handle_duplicate_features(features_df)
        nan_count = features_df.isna().sum().sum()
        if nan_count > 0: features_df.fillna(0, inplace=True)
        print(f"Finished basic processing for RNA file {filename}.")
        return features_df
    except Exception as e: print(f"  Error during basic RNA processing for '{filename}': {e}"); return None


# --- normalize_rna_data function ---
def normalize_rna_data(rna_features_orig_index):
    """Filters low counts, normalizes (CPM+Log2). Input indexed by original_rna_sample_id."""
    print("\n--- Filtering and Normalizing Combined RNA Data ---")
    if rna_features_orig_index is None or rna_features_orig_index.empty: return None
    if rna_features_orig_index.index.name != 'original_rna_sample_id': print("ERROR: Input must be indexed by 'original_rna_sample_id'."); return None
    rna_counts = rna_features_orig_index.select_dtypes(include=np.number).fillna(0).astype(np.float32)
    if rna_counts.empty: print("Error: No numeric gene columns."); return None
    print(f"Input shape: {rna_counts.shape}")
    min_samples = max(2, int(len(rna_counts) * RNA_FILTER_MIN_SAMPLES_FRAC))
    print(f"Filtering: Keeping genes >= {RNA_FILTER_MIN_COUNTS} counts in >= {min_samples} samples.")
    try: genes_to_keep_mask = (rna_counts.values >= RNA_FILTER_MIN_COUNTS).sum(axis=0) >= min_samples; rna_filtered = rna_counts.loc[:, genes_to_keep_mask]
    except Exception as e: print(f"Error during gene filtering mask: {e}"); return None
    if rna_filtered.empty: print("Error: Filtering removed all genes."); return None
    print(f"Removed {rna_counts.shape[1] - rna_filtered.shape[1]} genes. Shape after filtering: {rna_filtered.shape}")
    print("Normalizing (CPM + Log2)...")
    try:
        library_sizes_col = rna_filtered.values.sum(axis=1)[:, np.newaxis] + 1e-9
        cpm_array = (rna_filtered.values / library_sizes_col) * 1e6
        log2_cpm_array = np.log2(np.maximum(cpm_array, 0) + 1)
        rna_normalized_log2 = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
        print("Normalization complete. Final shape:", rna_normalized_log2.shape)
        return rna_normalized_log2
    except Exception as e: print(f"Error during normalization: {e}"); return None


# --- align_samples function (Expects original_rna_sample_id index) ---
def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
    """Aligns two DataFrames based on common indices (expects 'original_rna_sample_id')."""
    print(f"\n--- Aligning {df1_name} and {df2_name} to Common Samples ---")
    if df1 is None or df2 is None or df1.empty or df2.empty: print("Input DataFrame missing/empty."); return None, None
    if df1.index.name != 'original_rna_sample_id' or df2.index.name != 'original_rna_sample_id': print("Error: Inputs must be indexed by 'original_rna_sample_id'."); return None, None
    print(f"Input {df1_name} shape: {df1.shape}; Input {df2_name} shape: {df2.shape}")
    common_samples = df1.index.intersection(df2.index)
    print(f"Found {len(common_samples)} common samples based on original_rna_sample_id index.")
    if len(common_samples) == 0: print("CRITICAL WARNING: No common samples found!"); return None, None
    try:
        df1_matched = df1.loc[common_samples].copy()
        df2_aligned = df2.loc[common_samples].reindex(df1_matched.index).copy()
        if not df1_matched.index.equals(df2_aligned.index): raise ValueError("Final indices mismatch!")
        print(f"Aligned shapes: {df1_name} {df1_matched.shape}, {df2_name} {df2_aligned.shape}. Alignment successful.")
        return df1_matched, df2_aligned # Returns indexed by original_rna_sample_id
    except Exception as e: print(f"Error during sample alignment: {e}"); return None, None


# --- apply_variance_filter function ---
def apply_variance_filter(df, n_top_features, data_type_name="Data"):
    """Applies variance filtering to keep top N features."""
    # Input df will be indexed by original_rna_sample_id after alignment
    # Output df_filtered will also be indexed by original_rna_sample_id
    print(f"\n--- Applying Variance Filter to {data_type_name} ---")
    if df is None or df.empty: print("Input DataFrame empty."); return None
    print(f"Input shape: {df.shape}. Targeting top {n_top_features} features.")
    try:
        df_numeric = df.select_dtypes(include=np.number)
        if df_numeric.empty: print("Error: No numeric columns."); return None
        variances = df_numeric.var(axis=0).dropna()
        if variances.empty: print("Error: No variance calculated."); return None
        n_to_keep = min(n_top_features, len(variances))
        print(f"Selecting top {n_to_keep} features.")
        top_features = variances.nlargest(n_to_keep).index.tolist()
        df_filtered = df[top_features].copy()
        print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
        return df_filtered # Returns df indexed by original_rna_sample_id
    except Exception as e: print(f"Error during variance filtering: {e}"); return None

# --- REMOVED align_data_and_index function ---
# --- REMOVED _format_id_component and _create_composite_id helpers ---
# (As they are not needed with the original_rna_sample_id indexing strategy)