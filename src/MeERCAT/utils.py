# meercat/meercat/utils.py

import re
import pandas as pd
import numpy as np
import os
import time # Keep time if setup_paths or other future utils need it

def setup_paths(base_path, project_folder, spearman_subfolder, nmf_subfolder, plots_subfolder):
    """Creates directory structure relative to base_path and returns paths dictionary."""
    paths = {}
    # Base output directory for the project
    paths['base'] = os.path.abspath(os.path.join(base_path, project_folder))
    # Specific analysis output directories
    paths['spearman'] = os.path.join(paths['base'], spearman_subfolder)
    paths['nmf'] = os.path.join(paths['base'], nmf_subfolder)
    # Plot directories within analysis folders
    paths['spearman_plots'] = os.path.join(paths['spearman'], plots_subfolder)
    paths['nmf_plots'] = os.path.join(paths['nmf'], plots_subfolder)

    # Create all defined directories
    print(f"\n--- Ensuring Output Directories Exist ---")
    for key, path in paths.items():
        try:
            os.makedirs(path, exist_ok=True)
            # print(f"  Confirmed/Created: {path}") # Optional confirmation
        except OSError as e:
            print(f"WARNING: Could not create directory {path}: {e}")
            # Depending on severity, you might want to raise the error or handle it
    print(f"Base output directory: {paths['base']}")
    return paths

def identify_metabolite_columns(columns, prefix_patterns, suffix_regex):
    """Identifies metabolite columns based on prefixes OR suffixes."""
    metabolite_cols = []
    non_metabolite_cols = []
    prefix_patterns_lower = [p.lower() for p in prefix_patterns]
    for col in columns:
        if isinstance(col, str):
            col_lower = col.lower().strip()
            is_prefix_match = any(col_lower.startswith(p) for p in prefix_patterns_lower)
            is_suffix_match = bool(re.search(suffix_regex, col_lower))
            if is_prefix_match or is_suffix_match:
                metabolite_cols.append(col)
            else:
                non_metabolite_cols.append(col)
        else: # Handle non-string column names if necessary
             non_metabolite_cols.append(col)
    return metabolite_cols, non_metabolite_cols

def clean_feature_names(columns):
    """Cleans feature IDs (e.g., removes .version numbers from gene IDs)."""
    original_cols = list(columns)
    # Ensure conversion to string before applying string methods
    cleaned_cols = pd.Index(original_cols).astype(str).str.replace(r'\.\d+$', '', regex=True)
    num_renamed = sum(1 for orig, clean in zip(original_cols, cleaned_cols) if orig != clean)
    if num_renamed > 0:
        print(f"  Cleaned {num_renamed} feature column names (removed trailing .version).")
    return cleaned_cols.tolist()

def handle_duplicate_features(df):
    """Handles duplicate feature columns by summing."""
    # Ensure columns are strings before checking duplicates
    df.columns = df.columns.astype(str)
    duplicated_features = df.columns[df.columns.duplicated(keep=False)]
    unique_duplicates = duplicated_features.unique()
    if not unique_duplicates.empty:
        print(f"  Found {len(unique_duplicates)} duplicate feature columns. Consolidating by summing...")
        # Ensure numeric before summing
        numeric_cols = df.select_dtypes(include=np.number).columns
        df_numeric_part = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
        # Group by column names (level=0 on axis=1) and sum
        df_consolidated = df_numeric_part.groupby(level=0, axis=1).sum()
        # Only return the consolidated numeric part
        df = df_consolidated
        print(f"  Shape after consolidating duplicates: {df.shape}")
    return df

# clean_index function is NOT in this file; it's moved to load_data.py
# Data loading functions (load_metadata, load_metabolite_files, etc.) are NOT in this file; they are in load_data.py