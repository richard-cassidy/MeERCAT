# meercat_analyzer/meercat_analyzer/preprocess.py

import pandas as pd
import numpy as np
import re
from .utils import identify_metabolite_columns, clean_feature_names, handle_duplicate_features, clean_index
from .config import METADATA_SAMPLE_COL, SAMPLE_FILTER_EXCLUDE_PATTERNS, \
                   METADATA_COMPOSITE_ID_COLS, RNA_FILTER_MIN_COUNTS, \
                   RNA_FILTER_MIN_SAMPLES_FRAC, VAR_FILTER_TOP_N_GENES, \
                   VAR_FILTER_TOP_N_METABOLITES, METABOLITE_ID_PREFIXES, \
                   METABOLITE_ID_SUFFIX_REGEX, RNA_METADATA_COLS_EXPECTED

def combine_dataframes(df_dict, axis=0, join='outer'):
    """Combines a dictionary of DataFrames."""
    if not df_dict:
        print("Warning: No DataFrames provided for combination.")
        return None
    print(f"\n--- Combining {len(df_dict)} DataFrames (axis={axis}, join='{join}') ---")
    try:
        combined_df = pd.concat(df_dict.values(), axis=axis, join=join, ignore_index=(axis==0))
        print(f"Combined shape: {combined_df.shape}")
        return combined_df
    except Exception as e:
        print(f"Error during concatenation: {e}")
        return None

def clean_metabolite_data(combined_df):
    """Cleans combined metabolite data: numeric conversion, NaN handling, filtering, indexing."""
    print("\n--- Cleaning Combined Metabolite Data ---")
    if combined_df is None or combined_df.empty:
        print("Input DataFrame is None or empty. Skipping cleaning.")
        return None

    df_copy = combined_df.copy()

    # 1. Identify columns
    all_cols = df_copy.columns.tolist()
    metabolite_feature_cols, non_metabolite_cols = identify_metabolite_columns(
        all_cols, METABOLITE_ID_PREFIXES, METABOLITE_ID_SUFFIX_REGEX
    )
    print(f"Identified {len(metabolite_feature_cols)} potential metabolite feature columns.")
    print(f"Found {len(non_metabolite_cols)} assumed non-feature columns: {sorted([str(c) for c in non_metabolite_cols])}")

    # 2. Convert features to numeric, coerce errors
    print("Converting features to numeric...")
    converted_count = 0
    for col in metabolite_feature_cols:
        if col in df_copy.columns and not pd.api.types.is_numeric_dtype(df_copy[col]):
            original_dtype = df_copy[col].dtype
            df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce')
            if original_dtype != df_copy[col].dtype: converted_count += 1
    print(f"Numeric conversion attempted ({converted_count} columns potentially changed).")

    # 3. Drop all-NaN columns (among features AND others)
    print("Checking for all-NaN columns...")
    nan_only_cols = [col for col in df_copy.columns if df_copy[col].isnull().all()]
    if nan_only_cols:
        print(f"Identified {len(nan_only_cols)} all-NaN columns to remove: {nan_only_cols}")
        df_copy.drop(columns=nan_only_cols, inplace=True)
        print(f"Dropped all-NaN columns. Shape now: {df_copy.shape}")
    else:
        print("No all-NaN columns found.")

    # 4. Filter out QC/Blank samples
    print("Filtering out QC/Blank samples...")
    if METADATA_SAMPLE_COL in df_copy.columns:
        exclude_regex = '|'.join(SAMPLE_FILTER_EXCLUDE_PATTERNS)
        # Ensure column is string before using .str
        sample_col_str = df_copy[METADATA_SAMPLE_COL].astype(str)
        keep_mask = ~sample_col_str.str.lower().str.contains(exclude_regex, na=False, regex=True)
        rows_before = len(df_copy)
        df_filtered = df_copy[keep_mask].copy()
        print(f"Filtered {rows_before - len(df_filtered)} QC/Blank/Other rows. Shape now: {df_filtered.shape}")
    else:
        print(f"Warning: '{METADATA_SAMPLE_COL}' column not found for filtering.")
        df_filtered = df_copy # Keep all rows

    if df_filtered.empty:
        print("Warning: DataFrame empty after filtering samples.")
        return df_filtered # Return empty df

    # 5. Create and Set Composite Index
    print("Creating composite sample index...")
    index_set_successfully = False
    missing_cols = [col for col in METADATA_COMPOSITE_ID_COLS if col not in df_filtered.columns]
    if missing_cols:
        print(f"Error: Missing required columns for composite index: {missing_cols}")
    else:
        try:
            # Build index components safely
            id_parts = []
            for col in METADATA_COMPOSITE_ID_COLS:
                 # Convert to numeric if possible, then format, handle NA
                 num_series = pd.to_numeric(df_filtered[col], errors='coerce')
                 str_series = num_series.apply(lambda x: str(int(x)) if pd.notna(x) and x == int(x) else str(x) if pd.notna(x) else 'na')
                 id_parts.append(str_series)

            composite_id_series = pd.Series("_".join(map(str, tpl)) for tpl in zip(*id_parts))
            composite_id_series.index = df_filtered.index # Align index before assigning

            # Check uniqueness BEFORE assigning as index
            if composite_id_series.duplicated().any():
                print("CRITICAL WARNING: Duplicate composite IDs generated! Index cannot be set.")
                # Store as column for inspection
                df_filtered['composite_sample_id_TEMP'] = composite_id_series
            else:
                print("Composite ID is unique. Setting as index...")
                df_filtered.index = clean_index(composite_id_series)
                df_filtered.index.name = 'composite_sample_id'
                print(f"Successfully set index: '{df_filtered.index.name}'")
                index_set_successfully = True

        except Exception as e:
            print(f"Error creating or setting composite index: {e}")

    print("Metabolite cleaning finished.")
    return df_filtered


def process_rna_dataframe(df_raw, filename, rows_are_genes=True):
    """Orients, cleans, labels, and indexes a single raw RNA DataFrame."""
    print(f"\nProcessing RNA DataFrame from '{filename}'...")
    if df_raw is None or df_raw.empty: return None

    try:
        # Extract Experiment ID
        match_exp = re.match(r'(MC\d+)', filename, re.IGNORECASE) or re.search(r'(MC\d+)', filename, re.IGNORECASE)
        experiment_id = match_exp.group(1).upper() if match_exp else "UNKNOWN_EXPERIMENT"

        # Identify potential sample/gene columns (heuristic)
        df_reset = df_raw.reset_index()
        potential_id_col = df_reset.columns[0]
        potential_sample_cols = [col for col in df_reset.columns if col != potential_id_col] # Simplistic guess

        # Reconstruct with index and separate numeric part
        df_with_index = df_reset.set_index(potential_id_col)
        df_numeric_part = df_with_index[potential_sample_cols] # Assume these are counts

        # Orient: samples as rows, genes as columns
        df_oriented = df_numeric_part.T.copy() if rows_are_genes else df_numeric_part.copy()
        print(f"  Oriented shape (Samples x Features): {df_oriented.shape}")

        # Clean index (Sample IDs) and columns (Gene IDs)
        df_oriented.index = clean_index(df_oriented.index)
        df_oriented.columns = clean_index(df_oriented.columns) # Clean gene IDs too

        # Clean Gene IDs (e.g., remove .version) & Handle Duplicates
        df_oriented.columns = clean_feature_names(df_oriented.columns)
        df_oriented = handle_duplicate_features(df_oriented)

        # --- Extract Metadata from Original Sample Index ---
        # This part is highly specific to the index format like "1_00_7days"
        print("  Extracting metadata from original sample index (heuristic)...")
        original_index_series = pd.Series(df_oriented.index) # Use the cleaned index now
        # Adjust pattern if needed - this is complex and error-prone
        pattern = r"(\d+)_(\d{2,})_(\d+)([a-zA-Z]*)"
        extracted_data = original_index_series.str.extract(pattern, expand=True)
        if extracted_data.shape[1] == 4:
            extracted_data.columns = ['replicate_str', 'concentration_str', 'duration_num_str', 'duration_unit']
            replicates = pd.to_numeric(extracted_data['replicate_str'], errors='coerce').astype('Int64')
            # Concentration mapping/conversion (add more cases if needed)
            conc_map = {'00': 0.0, '05': 0.5, '20': 2.0, '40': 4.0}
            concentrations_numeric = extracted_data['concentration_str'].map(conc_map)
            needs_direct_conv = concentrations_numeric.isnull() & extracted_data['concentration_str'].notnull()
            if needs_direct_conv.any():
                 concentrations_numeric[needs_direct_conv] = pd.to_numeric(extracted_data.loc[needs_direct_conv, 'concentration_str'], errors='coerce')
            # Day conversion
            days_numeric = pd.to_numeric(extracted_data['duration_num_str'], errors='coerce')
            is_weeks = extracted_data['duration_unit'].str.lower().str.contains('week', na=False)
            if is_weeks.any(): days_numeric[is_weeks] = days_numeric[is_weeks] * 7
        else:
            print("  WARNING: Could not extract metadata reliably from index pattern. Metadata columns will be NaN.")
            replicates = pd.Series(np.nan, index=df_oriented.index)
            concentrations_numeric = pd.Series(np.nan, index=df_oriented.index)
            days_numeric = pd.Series(np.nan, index=df_oriented.index)

        # Create treatment group
        treatment_groups = np.where(concentrations_numeric == 0.0, "CONTROL_0.0",
            np.where(pd.notna(concentrations_numeric), "ARSENIC_" + concentrations_numeric.astype(str), "UNKNOWN_TREATMENT"))

        # --- Create Composite ID ---
        print("  Creating composite sample ID...")
        exp_id_str_rna = str(experiment_id); conc_str_rna = concentrations_numeric.astype(str)
        day_str_rna = pd.Series(days_numeric).apply(lambda x: str(int(x)) if pd.notna(x) and x==int(x) else str(x) if pd.notna(x) else 'na')
        rep_str_rna = replicates.astype(str)
        composite_id_series = (exp_id_str_rna + "_" + conc_str_rna.fillna('na') + "_" + day_str_rna.fillna('na') + "_" + rep_str_rna.fillna('na'))
        composite_id_col_name = 'composite_sample_id'
        df_oriented[composite_id_col_name] = composite_id_series.values

        # --- Set Composite Index ---
        if df_oriented[composite_id_col_name].duplicated().any():
             print(f"  WARNING: Duplicate composite IDs generated within RNA file {filename}! Index not set.")
             # Keep metadata as columns
             df_oriented['replicate'] = replicates.values; df_oriented['arsenic_concentration'] = concentrations_numeric.values
             df_oriented['days'] = days_numeric.values; df_oriented['treatment_group'] = treatment_groups; df_oriented['experiment_id'] = experiment_id
        else:
             df_oriented.set_index(composite_id_col_name, inplace=True)
             df_oriented.index = clean_index(df_oriented.index) # Clean final index
             print(f"  Set composite ID as index. Example: {df_oriented.index[0]}")
             # Add metadata as columns AFTER setting index
             df_oriented['replicate'] = replicates.values; df_oriented['arsenic_concentration'] = concentrations_numeric.values
             df_oriented['days'] = days_numeric.values; df_oriented['treatment_group'] = treatment_groups; df_oriented['experiment_id'] = experiment_id

        return df_oriented

    except Exception as e:
         print(f"  An unexpected error occurred during RNA processing for '{filename}': {e}")
         return None

def normalize_rna_data(combined_raw_counts_df):
    """Filters low-count genes and performs CPM + log2(x+1) normalization."""
    print("\n--- Filtering and Normalizing RNA Data ---")
    if combined_raw_counts_df is None or combined_raw_counts_df.empty:
        print("Input DataFrame is None or empty. Skipping.")
        return None

    # Separate metadata (assuming it exists as columns)
    metadata_cols = [col for col in combined_raw_counts_df.columns if col in RNA_METADATA_COLS_EXPECTED]
    rna_metadata = combined_raw_counts_df[metadata_cols].copy() if metadata_cols else pd.DataFrame(index=combined_raw_counts_df.index)
    gene_cols = [col for col in combined_raw_counts_df.columns if col not in metadata_cols]

    if not gene_cols:
        print("Error: No gene count columns identified.")
        return None

    rna_counts = combined_raw_counts_df[gene_cols].apply(pd.to_numeric, errors='coerce').fillna(0).astype(np.float32)
    print(f"Separated counts. Shape (Samples x Genes): {rna_counts.shape}")

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

    rna_normalized_log2_counts = pd.DataFrame(log2_cpm_array, index=rna_filtered.index, columns=rna_filtered.columns)
    print("Normalization complete.")

    # 3. Re-attach Metadata
    print("Re-attaching metadata...")
    # Align metadata index to normalized counts index
    rna_metadata_aligned = rna_metadata.reindex(rna_normalized_log2_counts.index)
    rna_normalized_log2_final = pd.concat([rna_metadata_aligned, rna_normalized_log2_counts], axis=1)

    if not rna_normalized_log2_final.index.equals(rna_metadata_aligned.index):
        print("CRITICAL Error: Index misalignment after re-attaching metadata.")
        return None
    else:
        print("Metadata successfully re-attached.")
        print("Final processed RNA data shape:", rna_normalized_log2_final.shape)
        return rna_normalized_log2_final


def align_samples(df1, df2, df1_name="DataFrame1", df2_name="DataFrame2"):
    """Aligns two DataFrames based on common indices."""
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
             return pd.DataFrame(index=df.index) # Return empty with index

        variances = df_numeric.var(axis=0).dropna() # Drop NaNs if variance is NaN
        if variances.empty:
             print("Error: Could not calculate variance for any feature.")
             return pd.DataFrame(index=df.index)

        n_to_keep = min(n_top_features, len(variances))
        print(f"Selecting top {n_to_keep} features (out of {len(variances)} with calculated variance).")
        top_features = variances.nlargest(n_to_keep).index.tolist()

        # Filter the original DataFrame (including potential non-numeric ID cols if they were kept)
        # using the selected feature names
        df_filtered = df[top_features].copy()
        print(f"Filtered {data_type_name} shape: {df_filtered.shape}")
        return df_filtered

    except Exception as e:
        print(f"Error during variance filtering for {data_type_name}: {e}")
        return None