# mercat_analyzer/meercat/preprocess.py

import pandas as pd
import numpy as np
import re
# Import utils and config functions/variables
from .utils import identify_metabolite_columns, clean_feature_names, handle_duplicate_features # REMOVED clean_index
from .config import (
    METADATA_SAMPLE_COL, SAMPLE_FILTER_EXCLUDE_PATTERNS,
    METADATA_COMPOSITE_ID_COLS, RNA_FILTER_MIN_COUNTS,
    RNA_FILTER_MIN_SAMPLES_FRAC, VAR_FILTER_TOP_N_GENES,
    VAR_FILTER_TOP_N_METABOLITES, METABOLITE_ID_PREFIXES,
    METABOLITE_ID_SUFFIX_REGEX, RNA_METADATA_COLS_EXPECTED
)
# **** ADDED: Import clean_index from load_data ****
try:
    from .load_data import clean_index
except ImportError:
    # Fallback or raise error if load_data structure is wrong
    print("ERROR in preprocess.py: Could not import clean_index from .load_data.")
    # Define a dummy function to avoid crashing later if needed for testing,
    # but this indicates a setup problem.
    def clean_index(index):
        print("ERROR: clean_index function missing!")
        return index # Return original index

# --- Helper for Consistent ID String Formatting ---
def _format_id_component(series):
    """Converts a numeric series to string for IDs, handling NaN and '.0'."""
    num_series = pd.to_numeric(series, errors='coerce')
    str_series = num_series.apply(
        lambda x: str(int(x)) if pd.notna(x) and x == int(x) else str(x) if pd.notna(x) else 'na'
    )
    return str_series

# --- Helper to Create Composite ID ---
def _create_composite_id(df, id_cols, exp_id=None):
    """Creates a composite ID Series from specified columns in a DataFrame."""
    print(f"  Attempting to create composite ID using columns: {id_cols}")
    id_parts = []
    # --- Handle Experiment ID ---
    if exp_id is not None:
         id_parts.append(pd.Series(str(exp_id), index=df.index)); print(f"    Using provided experiment_id: {exp_id}")
         id_cols = [col for col in id_cols if col != 'experiment_id'] # Avoid double processing
    elif 'experiment_id' in id_cols and 'experiment_id' in df.columns:
         id_parts.append(df['experiment_id'].astype(str).fillna('na')); print("    Using 'experiment_id' column.")
         id_cols = [col for col in id_cols if col != 'experiment_id']
    else: print("  ERROR: Experiment ID missing."); return None
    # --- Check remaining ---
    missing_cols = [col for col in id_cols if col not in df.columns]
    if missing_cols: print(f"  ERROR: Missing columns for composite ID: {missing_cols}"); return None
    # --- Format other components ---
    print(f"    Processing components: {id_cols}")
    for col in id_cols: id_parts.append(_format_id_component(df[col]))
    # --- Combine ---
    try:
        if not id_parts: print("  ERROR: No components for composite ID."); return None
        composite_id_series = pd.Series("_".join(map(str, tpl)) for tpl in zip(*id_parts))
        composite_id_series.index = df.index; print("  Composite ID Series created successfully.")
        return composite_id_series
    except Exception as e_zip: print(f"  ERROR during ID component combination: {e_zip}"); return None

# --- combine_dataframes function ---
def combine_dataframes(df_dict, axis=0, join='outer'):
    """Combines a dictionary of DataFrames."""
    if not df_dict: print("Warning: No DataFrames provided for combination."); return None
    print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
    try:
        valid_dfs = [df for df in df_dict.values() if isinstance(df, pd.DataFrame)]
        if not valid_dfs: print("Warning: No valid DataFrames in dictionary values."); return None
        combined_df = pd.concat(valid_dfs, axis=axis, join=join, ignore_index=(axis==0))
        print(f"Combined shape: {combined_df.shape}")
        return combined_df
    except Exception as e: print(f"Error during concatenation: {e}"); return None

# --- clean_metabolite_data function ---
def clean_metabolite_data(combined_df):
    """Cleans combined metabolite data: separates metadata, cleans features, filters, creates index."""
    print("\n--- Cleaning Combined Metabolite Data ---")
    if combined_df is None or combined_df.empty: return None, None
    df_copy = combined_df.copy(); metadata_df = None; index_set_successfully = False
    all_cols = df_copy.columns.tolist()
    metabolite_feature_cols, non_metabolite_cols = identify_metabolite_columns(all_cols, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX)
    print(f"Identified {len(metabolite_feature_cols)} potential features.")
    potential_metadata_cols = list(set(non_metabolite_cols + METADATA_COMPOSITE_ID_COLS + [METADATA_SAMPLE_COL]))
    metadata_cols_present = [col for col in potential_metadata_cols if col in df_copy.columns]
    if metadata_cols_present: metadata_df = df_copy[metadata_cols_present].copy()
    else: metadata_df = pd.DataFrame(index=df_copy.index)
    print("Converting features to numeric...")
    feature_cols_present = [col for col in metabolite_feature_cols if col in df_copy.columns]
    df_features = df_copy[feature_cols_present].copy()
    for col in feature_cols_present:
        if not pd.api.types.is_numeric_dtype(df_features[col]): df_features[col] = pd.to_numeric(df_features[col], errors='coerce')
    numeric_feature_cols = df_features.select_dtypes(include=np.number).columns.tolist()
    df_numeric_features = df_features[numeric_feature_cols].copy()
    print(f"Kept {len(numeric_feature_cols)} numeric feature columns.")
    nan_feature_cols = [col for col in numeric_feature_cols if df_numeric_features[col].isnull().all()]
    if nan_feature_cols: df_numeric_features.drop(columns=nan_feature_cols, inplace=True); print(f"Dropped {len(nan_feature_cols)} all-NaN feature columns.")
    df_to_filter = pd.merge(metadata_df, df_numeric_features, left_index=True, right_index=True, how='inner')
    print("Filtering out QC/Blank samples...")
    df_filtered = df_to_filter
    if METADATA_SAMPLE_COL in df_filtered.columns:
        exclude_regex = '|'.join(SAMPLE_FILTER_EXCLUDE_PATTERNS)
        keep_mask = ~df_filtered[METADATA_SAMPLE_COL].astype(str).str.lower().str.contains(exclude_regex, na=False, regex=True)
        df_filtered = df_filtered[keep_mask].copy()
        print(f"Shape after filtering samples: {df_filtered.shape}")
    if df_filtered.empty: print("Warning: DataFrame empty after filtering."); return df_filtered, pd.DataFrame()
    final_metadata_cols = [col for col in metadata_cols_present if col in df_filtered.columns]
    metadata_df_final = df_filtered[final_metadata_cols].copy()
    final_feature_cols = [col for col in df_filtered.columns if col not in final_metadata_cols]
    df_features_only = df_filtered[final_feature_cols].copy()
    print("Creating composite sample index...")
    composite_id_series = _create_composite_id(metadata_df_final, METADATA_COMPOSITE_ID_COLS)
    if composite_id_series is not None:
        if composite_id_series.duplicated().any():
            print("CRITICAL WARNING: Duplicate metabolite composite IDs! Index not set."); metadata_df_final['composite_id_TEMP'] = composite_id_series; df_features_only['composite_id_TEMP'] = composite_id_series
        else:
            print("Setting metabolite composite index...")
            cleaned_final_index = clean_index(composite_id_series) # Uses clean_index
            try: metadata_df_final.index = cleaned_final_index; metadata_df_final.index.name = 'composite_sample_id'; df_features_only.index = cleaned_final_index; df_features_only.index.name = 'composite_sample_id'; index_set_successfully = True; print("Index set successfully.")
            except Exception as e_idx: print(f"ERROR setting index: {e_idx}")
    else: print("Failed to create metabolite composite ID.")
    print("Metabolite cleaning finished.")
    return df_features_only, metadata_df_final # Return regardless of index success for now

# --- process_rna_dataframe function ---
def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
    """Orients, cleans features, returns RNA features indexed by original sample ID."""
    print(f"\nProcessing RNA DataFrame from '{filename}'...")
    if df_raw is None or df_raw.empty: return None
    try:
        if rows_are_genes:
            df_reset = df_raw.reset_index(); gene_id_col = df_reset.columns[0]
            sample_cols = [col for col in df_reset.columns if col != gene_id_col]
            if not sample_cols: raise ValueError("No sample columns found.")
            df_numeric_part = df_reset.set_index(gene_id_col)[sample_cols].apply(pd.to_numeric, errors='coerce')
            features_df = df_numeric_part.T.copy()
        else: features_df = df_raw.apply(pd.to_numeric, errors='coerce').copy()
        print(f"  Oriented shape (Samples x Features): {features_df.shape}")
        features_df.index = clean_index(features_df.index) # Uses clean_index
        features_df.index.name = 'original_rna_sample_id'
        features_df.columns = clean_index(features_df.columns) # Uses clean_index
        features_df.columns = clean_feature_names(features_df.columns)
        features_df = handle_duplicate_features(features_df)
        nan_count = features_df.isna().sum().sum()
        if nan_count > 0: features_df.fillna(0, inplace=True); print(f"  Filled {nan_count} NaNs with 0.")
        print(f"Finished basic processing for RNA file {filename}.")
        return features_df
    except Exception as e: print(f"  Error during basic RNA processing for '{filename}': {e}"); return None

# --- normalize_rna_data function ---
def normalize_rna_data(rna_features_indexed):
    """Filters low-count genes and performs CPM + log2(x+1) normalization."""
    print("\n--- Filtering and Normalizing Combined RNA Data ---")
    if rna_features_indexed is None or rna_features_indexed.empty: return None
    rna_counts = rna_features_indexed.select_dtypes(include=np.number).fillna(0).astype(np.float32)
    if rna_counts.empty: print("Error: No numeric gene columns for normalization."); return None
    print(f"Input shape: {rna_counts.shape}")
    min_samples = max(2, int(len(rna_counts) * RNA_FILTER_MIN_SAMPLES_FRAC))
    print(f"Filtering: Keeping genes >= {RNA_FILTER_MIN_COUNTS} counts in >= {min_samples} samples.")
    try:
        genes_to_keep_mask = (rna_counts.values >= RNA_FILTER_MIN_COUNTS).sum(axis=0) >= min_samples
        rna_filtered = rna_counts.loc[:, genes_to_keep_mask]
        print(f"Removed {rna_counts.shape[1] - rna_filtered.shape[1]} genes. Shape after filtering: {rna_filtered.shape}")
    except Exception as e: print(f"Error during gene filtering mask: {e}"); return None
    if rna_filtered.empty: print("Error: Filtering removed all genes."); return None
    print("Normalizing (CPM + Log2)...")
    try:
        library_sizes_col = rna_filtered.values.sum(axis=1)[:, np.newaxis] + 1e-9
        cpm_array = (rna_filtered.values / library_sizes_col) * 1e6
        log2_cpm_array = np.log2(np.maximum(cpm_array, 0) + 1)
        rna_normalized_log2 = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
        print("Normalization complete. Final shape:", rna_normalized_log2.shape)
        return rna_normalized_log2
    except Exception as e: print(f"Error during normalization: {e}"); return None

# --- align_samples function ---
def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
    """Aligns two DataFrames based on common indices (assumes 'composite_sample_id')."""
    print(f"\n--- Aligning {df1_name} and {df2_name} to Common Samples ---")
    if df1 is None or df2 is None or df1.empty or df2.empty: print("Input DataFrame missing/empty."); return None, None
    if df1.index.name != 'composite_sample_id' or df2.index.name != 'composite_sample_id': print("Error: Input DFs must be indexed by 'composite_sample_id'."); return None, None
    print(f"Input {df1_name} shape: {df1.shape}; Input {df2_name} shape: {df2.shape}")
    # Indices should be clean already, but apply clean_index again for safety
    index1_clean = clean_index(df1.index) # Uses clean_index
    index2_clean = clean_index(df2.index) # Uses clean_index
    common_samples = index1_clean.intersection(index2_clean)
    print(f"Found {len(common_samples)} common samples.")
    if len(common_samples) == 0: print("CRITICAL WARNING: No common samples found!"); return None, None
    try:
        df1_matched = df1.loc[common_samples].copy(); df2_matched = df2.loc[common_samples].copy()
        df1_matched.index = clean_index(df1_matched.index); df2_matched.index = clean_index(df2_matched.index) # Re-clean after loc
        df2_aligned = df2_matched.reindex(df1_matched.index)
        if not df1_matched.index.equals(df2_aligned.index): raise ValueError("Final indices mismatch!")
        print(f"Aligned shapes: {df1_name} {df1_matched.shape}, {df2_name} {df2_aligned.shape}. Alignment successful.")
        return df1_matched, df2_aligned
    except Exception as e: print(f"Error during sample alignment: {e}"); return None, None

# --- apply_variance_filter function ---
def apply_variance_filter(df, n_top_features, data_type_name="Data"):
    """Applies variance filtering to keep top N features."""
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
        df_filtered = df[top_features].copy() # Keep only selected features
        print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
        return df_filtered
    except Exception as e: print(f"Error during variance filtering: {e}"); return None