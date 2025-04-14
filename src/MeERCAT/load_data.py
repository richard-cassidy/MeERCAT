# mercat_analyzer/meercat/load_data.py

import pandas as pd
import os
import io
import re
# NOTE: No relative import needed as clean_index is defined below

# ==================================================
# Helper Function(s) needed within this module
# ==================================================
def clean_index(index):
    """Converts index to string, lowercases, and strips whitespace."""
    # Check if input is already a pandas Index
    if isinstance(index, pd.Index):
        # Ensure it's string type before using .str methods
        return index.astype(str).str.lower().str.strip()
    # Handle other potential iterable inputs (like lists, Series)
    elif isinstance(index, (list, pd.Series)):
         # print("Warning: Input is not a pandas Index, converting.") # Less verbose
         return pd.Index(index).astype(str).str.lower().str.strip()
    # Handle single string or other types gracefully
    elif isinstance(index, str):
         return index.lower().strip()
    else:
         # print(f"Warning: Unexpected input type ({type(index)}) for clean_index.") # Less verbose
         try:
             # Attempt conversion assuming input might be a single item not in a list
             return pd.Index([str(index)]).astype(str).str.lower().str.strip()[0]
         except Exception as e:
             print(f"ERROR: Could not convert input to cleanable index: {e}")
             return index # Return original if conversion fails

# ==================================================
# Data Loading Functions
# ==================================================

def load_metadata(metadata_filepath):
    """Loads metadata from a CSV file."""
    print(f"\n--- Loading Metadata from: {metadata_filepath} ---")
    if not os.path.exists(metadata_filepath):
        print(f"ERROR: Metadata file not found at '{metadata_filepath}'")
        return None
    try:
        metadata_df = pd.read_csv(metadata_filepath)
        print(f"Metadata loaded successfully. Shape: {metadata_df.shape}")
        # print("Head:\n", metadata_df.head()) # Optional
        return metadata_df
    except Exception as e:
        print(f"Error loading metadata file: {e}")
        return None

# --- Load Metabolites ---
def load_metabolite_files(metabolite_directory_path):
    """Loads multiple metabolite CSV files from a directory, standardizes columns, adds IDs."""
    print(f"\n--- Loading Metabolite Files from Directory: {metabolite_directory_path} ---")
    processed_data = {} # Dictionary to store loaded DataFrames keyed by filename
    if not os.path.isdir(metabolite_directory_path):
        print(f"ERROR: Metabolite input directory not found: {metabolite_directory_path}")
        return processed_data # Return empty dict
    metabolite_file_dict = { f: os.path.join(metabolite_directory_path, f)
                             for f in os.listdir(metabolite_directory_path) if f.lower().endswith('.csv') }
    if not metabolite_file_dict:
        print(f"Warning: No CSV files found in directory: {metabolite_directory_path}")
        return processed_data # Return empty dict
    print(f"Found {len(metabolite_file_dict)} potential CSV files to load.")

    for filename, filepath in metabolite_file_dict.items():
        print(f"\nProcessing: {filename}")
        try:
            df = pd.read_csv(filepath)
            df.columns = [str(c).lower().strip() for c in df.columns] # Standardize cols
            match = re.match(r'(MC\d+)', filename, re.IGNORECASE) or re.search(r'(MC\d+)', filename, re.IGNORECASE)
            experiment_id = match.group(1).upper() if match else "UNKNOWN_EXPERIMENT"
            if experiment_id == "UNKNOWN_EXPERIMENT": print(f"  Warning: No MC ID found in filename '{filename}'.")
            df['experiment_id'] = experiment_id
            df['source_filename'] = filename
            processed_data[filename] = df
        except Exception as e: print(f"  Error processing file {filename}: {e}")

    if not processed_data: print("\nERROR: No metabolite files were successfully loaded.")
    else: print(f"\nFinished loading metabolite files. {len(processed_data)} DataFrames loaded.")
    return processed_data


# --- Load RNA ---
def load_rna_files(rna_directory_path, rows_are_genes=True):
    """Loads multiple RNA seq count CSV files from a directory."""
    print(f"\n--- Loading RNA Count Files from Directory: {rna_directory_path} ---")
    loaded_data = {}
    if not os.path.isdir(rna_directory_path):
        print(f"ERROR: RNA input directory not found: {rna_directory_path}")
        return loaded_data
    rna_file_dict = { f: os.path.join(rna_directory_path, f)
                      for f in os.listdir(rna_directory_path) if f.lower().endswith('.csv')}
    if not rna_file_dict:
        print(f"Warning: No CSV files found in directory: {rna_directory_path}")
        return loaded_data
    print(f"Found {len(rna_file_dict)} potential CSV files to load.")

    for filename, filepath in rna_file_dict.items():
        print(f"\nProcessing RNA file: {filename}")
        try:
            # Load as string initially, specify index column
            # Use sep=None for auto-detection, engine='python' for flexibility
            df_raw = pd.read_csv(filepath, index_col=0, sep=None, engine='python', dtype=str)
            print(f"  Loaded '{filename}'. Initial shape: {df_raw.shape}")
            loaded_data[filename] = df_raw
        except Exception as e: print(f"  Error loading RNA count file '{filename}': {e}")

    if not loaded_data: print("\nERROR: No RNA files were successfully loaded.")
    else: print(f"\nFinished loading RNA files. {len(loaded_data)} DataFrames loaded.")
    return loaded_data


# --- Load RNA Metadata (Uses the local clean_index) ---
def load_external_rna_metadata(rna_metadata_filepath):
    """Loads external RNA metadata CSV, indexed by cleaned original_rna_sample_id."""
    print(f"\n--- Loading External RNA Metadata: {os.path.basename(rna_metadata_filepath)} ---")

    # 1. Check file existence
    if not os.path.exists(rna_metadata_filepath):
        print(f"ERROR: File not found at '{rna_metadata_filepath}'")
        return None

    # 2. Try loading, indexing, and cleaning
    try:
        # Assume first col is index ('original_rna_sample_id')
        metadata_df = pd.read_csv(rna_metadata_filepath, index_col=0)

        # Check if index was actually created
        if not isinstance(metadata_df.index, pd.Index):
             raise TypeError("Loaded data does not have a valid index after read_csv (check index_col=0).")

        # *** Call the clean_index function defined WITHIN this file ***
        metadata_df.index = clean_index(metadata_df.index)
        metadata_df.index.name = 'original_rna_sample_id' # Ensure index name is set

        print(f"  Successfully loaded and indexed. Shape: {metadata_df.shape}")
        return metadata_df

    except Exception as e: # Catch general exceptions during loading/indexing
        print(f"ERROR loading or processing RNA metadata file '{os.path.basename(rna_metadata_filepath)}': {e}")
        return None