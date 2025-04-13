# mercat_analyzer/mercat_analyzer/load_data.py

import pandas as pd
import os
import io
import re

def load_metadata(metadata_filepath):
    """Loads metadata from a CSV file."""
    print(f"\n--- Loading Metadata from: {metadata_filepath} ---")
    if not os.path.exists(metadata_filepath):
        print(f"ERROR: Metadata file not found at '{metadata_filepath}'")
        return None
    try:
        metadata_df = pd.read_csv(metadata_filepath)
        print(f"Metadata loaded successfully. Shape: {metadata_df.shape}")
        print("Head:\n", metadata_df.head())
        return metadata_df
    except Exception as e:
        print(f"Error loading metadata file: {e}")
        return None

def load_metabolite_files(metabolite_file_dict):
    """Loads multiple metabolite CSV files, standardizes columns, adds IDs."""
    print(f"\n--- Loading {len(metabolite_file_dict)} Metabolite Files ---")
    processed_data = {}
    if not metabolite_file_dict:
        print("Warning: No metabolite files provided in dictionary.")
        return processed_data

    for filename, filepath in metabolite_file_dict.items():
        print(f"\nProcessing: {filename} ({filepath})")
        if not os.path.exists(filepath):
            print("  ERROR: File not found.")
            continue
        if not filename.lower().endswith('.csv'):
            print("  ERROR: File is not a CSV. Skipping.")
            continue

        try:
            # Load data
            df = pd.read_csv(filepath)
            print(f"  Loaded. Initial shape: {df.shape}")

            # Standardize column names
            df.columns = [str(c).lower().strip() for c in df.columns]

            # Extract Experiment ID
            match = re.match(r'(MC\d+)', filename, re.IGNORECASE) or re.search(r'(MC\d+)', filename, re.IGNORECASE)
            experiment_id = match.group(1).upper() if match else "UNKNOWN_EXPERIMENT"
            if experiment_id == "UNKNOWN_EXPERIMENT": print(f"  Warning: No MC ID found in filename '{filename}'.")

            # Add ID columns
            df['experiment_id'] = experiment_id
            df['source_filename'] = filename
            processed_data[filename] = df
            print(f"  Added IDs and stored DataFrame.")

        except Exception as e:
            print(f"  Error processing file {filename}: {e}")

    print(f"\nFinished loading metabolite files. {len(processed_data)} DataFrames loaded.")
    return processed_data


def load_rna_files(rna_file_dict, rows_are_genes=True):
    """Loads multiple RNA seq count files (expects gene x sample or sample x gene)."""
    print(f"\n--- Loading {len(rna_file_dict)} RNA Count Files ---")
    loaded_data = {}
    if not rna_file_dict:
        print("Warning: No RNA files provided.")
        return loaded_data

    for filename, filepath in rna_file_dict.items():
        print(f"\nProcessing RNA file: {filename} ({filepath})")
        if not os.path.exists(filepath):
            print("  ERROR: File not found.")
            continue
        if not filename.lower().endswith('.csv'): # Basic check, could allow tsv etc.
            print("  ERROR: File is not a CSV. Skipping.")
            continue
        try:
            # Load as string initially to prevent issues, specify index column
            df_raw = pd.read_csv(filepath, index_col=0, sep=None, engine='python', dtype=str)
            print(f"  Loaded '{filename}'. Initial shape: {df_raw.shape}")
            loaded_data[filename] = df_raw
        except Exception as e:
            print(f"  Error loading RNA count file '{filename}': {e}")

    print(f"\nFinished loading RNA files. {len(loaded_data)} DataFrames loaded.")
    return loaded_data