# mercat_analyzer/mercat_analyzer/utils.py

import re
import pandas as pd
import numpy as np
import os
import time

def setup_paths(base_path, project_folder, spearman_subfolder, nmf_subfolder, plots_subfolder):
    """Creates directory structure and returns paths."""
    paths = {}
    paths['base'] = os.path.join(base_path, project_folder)
    paths['spearman'] = os.path.join(paths['base'], spearman_subfolder)
    paths['nmf'] = os.path.join(paths['base'], nmf_subfolder)
    paths['spearman_plots'] = os.path.join(paths['spearman'], plots_subfolder)
    paths['nmf_plots'] = os.path.join(paths['nmf'], plots_subfolder)

    for path in paths.values():
        try:
            os.makedirs(path, exist_ok=True)
        except OSError as e:
            print(f"Warning: Could not create directory {path}: {e}")
            # Depending on severity, you might want to raise the error
    print(f"Ensured directories exist under: {paths['base']}")
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
    """Cleans gene IDs (e.g., removes .version numbers)."""
    original_cols = list(columns)
    cleaned_cols = pd.Index(original_cols).astype(str).str.replace(r'\.\d+$', '', regex=True)
    num_renamed = sum(1 for orig, clean in zip(original_cols, cleaned_cols) if orig != clean)
    if num_renamed > 0:
        print(f"  Cleaned {num_renamed} feature column names (removed trailing .version).")
    return cleaned_cols.tolist()

def handle_duplicate_features(df):
    """Handles duplicate feature columns by summing."""
    duplicated_features = df.columns[df.columns.duplicated(keep=False)]
    unique_duplicates = duplicated_features.unique()
    if not unique_duplicates.empty:
        print(f"  Found {len(unique_duplicates)} duplicate feature columns. Consolidating by summing...")
        # Ensure numeric before summing
        numeric_cols = df.select_dtypes(include=np.number).columns
        non_numeric_cols = df.select_dtypes(exclude=np.number).columns
        df_numeric_part = df[numeric_cols].apply(pd.to_numeric, errors='coerce').fillna(0)
        df_consolidated = df_numeric_part.groupby(level=0, axis=1).sum()
        # Re-attach non-numeric if needed, though usually not needed for feature matrix
        # df = pd.concat([df[non_numeric_cols], df_consolidated], axis=1)
        df = df_consolidated # Usually just keep the numeric part
        print(f"  Shape after consolidating duplicates: {df.shape}")
    return df

def clean_index(index):
    """Converts index to string, lowercases, and strips whitespace."""
    if isinstance(index, pd.Index):
        return index.astype(str).str.lower().str.strip()
    else:
        print("Warning: Input is not a pandas Index. Attempting conversion.")
        return pd.Index(index).astype(str).str.lower().str.strip()