# Load Data Function
import pandas as pd
import os
from typing import Dict, Optional
from meercat import config
import pandas as pd
import os
from typing import Dict, Optional




# Helper Function 1: Load RNA and Metabolite Data sets
def load_data(input_dir: str) -> Dict[str, Dict[str, pd.DataFrame]]:
    rna_data = {}
    metabolite_data = {}
    if not os.path.isdir(input_dir):
        return {"rna_data": rna_data, "metabolite_data": metabolite_data}
    try:
        for filename in os.listdir(input_dir):
            file_path = os.path.join(input_dir, filename)
            if os.path.isfile(file_path):
                if filename.endswith('_RNA.csv'):
                    exp_id = filename.removesuffix('_RNA.csv')
                    if exp_id:
                        try: rna_data[exp_id] = pd.read_csv(file_path)
                        except Exception as e: print(f"  * Warning: Failed to load RNA file {filename}: {e}")
                elif filename.endswith('_metabolites.csv'):
                    exp_id = filename.removesuffix('_metabolites.csv')
                    if exp_id:
                       try: metabolite_data[exp_id] = pd.read_csv(file_path)
                       except Exception as e: print(f"  * Warning: Failed to load Metabolite file {filename}: {e}")
        print(f"   loaded RNA and metabolite files")
    except FileNotFoundError:
         print(f"  * Error: Input directory for omics data not found or became inaccessible: {input_dir}")
    except Exception as e:
        print(f"  * Error: An unexpected error occurred while reading omics data files: {e}")
    return {"rna_data": rna_data, "metabolite_data": metabolite_data}

# Helper Function 2: Load Metadata 
def load_metadata(
    input_dir: str,
    experiment_metadata_filename: str = config.EXPERIMENT_MATADATA_FILENAME,
    rna_metadata_filename: str = config.RNA_METADATA_FILENAME
) -> Dict[str, Optional[pd.DataFrame]]:
    metadata = {"experiment": None, "rna": None}
    if not os.path.isdir(input_dir):
         # Removed print here, handled by caller or main function check
         return metadata
    # Load experiment metadata
    try:
        exp_meta_path = os.path.join(input_dir, experiment_metadata_filename)
        metadata['experiment'] = pd.read_csv(exp_meta_path)
        print(f"   loaded experiment metadata")

    except FileNotFoundError:
         print(f"  * Warning: Experiment metadata file not found: {exp_meta_path}")
    except Exception as e:
        print(f"  * Warning: Failed to load Experiment metadata file {experiment_metadata_filename}: {e}")

    # Load RNA metadata
    try:
        rna_meta_path = os.path.join(input_dir, rna_metadata_filename)
        metadata['rna'] = pd.read_csv(rna_meta_path)
        print(f"   loaded RNA metadata")
    except FileNotFoundError:
        print(f"  * Warning: RNA metadata file not found: {rna_meta_path}")
    except Exception as e:
        print(f"  * Warning: Failed to load RNA metadata file {rna_metadata_filename}: {e}")
    return metadata


### Main Loading Function  ###
def load_all_data(
    input_dir: str,
    experiment_metadata_filename: str = config.EXPERIMENT_MATADATA_FILENAME,
    rna_metadata_filename: str = config.RNA_METADATA_FILENAME
) -> Dict:
    """
    Loads all data (metadata, RNA, metabolites) using helpers, adds
    simple warnings if specific data types were not loaded, and provides
    a summary list of successfully loaded data.
    
    Requirements:
        Input dir:
            - exeriment metadata 
            - rna metadata
            - metabolite data
            - rna data
    Args:
        input_dir (str): Path to the directory containing all input CSV files.
        experiment_metadata_filename (str): Filename for experiment metadata.
        rna_metadata_filename (str): Filename for RNA metadata.

    Returns:
        Dict: A dictionary containing:
            'experiment_metadata': DataFrame or None.
            'rna_metadata': DataFrame or None.
            'rna_data': Dict mapping experiment_id -> RNA DataFrame.
            'metabolite_data': Dict mapping experiment_id -> Metabolite DataFrame.
    """
    print ("\n\n\nLoading data...")
    #print(f"From directory: {input_dir}")

    # Check if input directory exists *before* calling helpers
    if not os.path.isdir(input_dir):
        print(f"CRITICAL ERROR: Input directory not found: {input_dir}")
        # Return an empty structure to avoid further errors
        return {
            "experiment_metadata": None, "rna_metadata": None,
            "rna_data": {}, "metabolite_data": {}
        }

    # Initialize the final dictionary structure
    all_data = {
        "experiment_metadata": None,
        "rna_metadata": None,
        "rna_data": {},
        "metabolite_data": {}
    }
    successfully_loaded_list = [] # List to store names of loaded items

    # --- Load Metadata ---
    #print("\nLoading metadata...")
    metadata_loaded = load_metadata(
        input_dir, experiment_metadata_filename, rna_metadata_filename
    )
    all_data["experiment_metadata"] = metadata_loaded.get('experiment')
    all_data["rna_metadata"] = metadata_loaded.get('rna')
    # Check results and update summary list
    if all_data["experiment_metadata"] is not None:
        successfully_loaded_list.append("Experiment Metadata")
    if all_data["rna_metadata"] is not None:
        successfully_loaded_list.append("RNA Metadata")
  


    # --- Load RNA/Metabolite Data ---
    #print("\nLoading RNA and Metabolite data...")
    omics_data_loaded = load_data(input_dir)
    all_data["rna_data"] = omics_data_loaded.get("rna_data", {})
    all_data["metabolite_data"] = omics_data_loaded.get("metabolite_data", {})
    # Check results and update summary list
    if all_data["rna_data"]: # Checks if the dictionary is not empty
        count = len(all_data['rna_data'])
        successfully_loaded_list.append(f"RNA Data ({count} experiments)")
    if all_data["metabolite_data"]: # Checks if the dictionary is not empty
        count = len(all_data['metabolite_data'])
        successfully_loaded_list.append(f"Metabolite Data ({count} experiments)")
    # else: # Optional warning
    #     print("  * Warning: No Metabolite data files were successfully loaded.")




    # --- Final Summary ---
    print("\nSuccessfully Loaded")
    if successfully_loaded_list:
        for item in successfully_loaded_list:
            print(f"  - {item}")
    else:
        print("  * Warning: No data components were successfully loaded.")
    print("----------------------------\n\n")


    return all_data


# # ==================================================
# # Helper Function(s) needed within this module
# # ==================================================
# def clean_index(index):
#     """Converts index to string, lowercases, and strips whitespace."""
#     # Check if input is already a pandas Index
#     if isinstance(index, pd.Index):
#         # Ensure it's string type before using .str methods
#         return index.astype(str).str.lower().str.strip()
#     # Handle other potential iterable inputs (like lists, Series)
#     elif isinstance(index, (list, pd.Series)):
#          # print("Warning: Input is not a pandas Index, converting.") # Less verbose
#          return pd.Index(index).astype(str).str.lower().str.strip()
#     # Handle single string or other types gracefully
#     elif isinstance(index, str):
#          return index.lower().strip()
#     else:
#          # print(f"Warning: Unexpected input type ({type(index)}) for clean_index.") # Less verbose
#          try:
#              # Attempt conversion assuming input might be a single item not in a list
#              return pd.Index([str(index)]).astype(str).str.lower().str.strip()[0]
#          except Exception as e:
#              print(f"ERROR: Could not convert input to cleanable index: {e}")
#              return index # Return original if conversion fails

# # ==================================================
# # Data Loading Functions
# # ==================================================

# def load_metadata(metadata_filepath):
#     """Loads metadata from a CSV file."""
#     print(f"\n--- Loading Metadata from: {metadata_filepath} ---")
#     if not os.path.exists(metadata_filepath):
#         print(f"ERROR: Metadata file not found at '{metadata_filepath}'")
#         return None
#     try:
#         metadata_df = pd.read_csv(metadata_filepath)
#         print(f"Metadata loaded successfully. Shape: {metadata_df.shape}")
#         # print("Head:\n", metadata_df.head()) # Optional
#         return metadata_df
#     except Exception as e:
#         print(f"Error loading metadata file: {e}")
#         return None





# # --- Load Metabolites ---
# def load_metabolite_files(metabolite_directory_path):
#     """Loads multiple metabolite CSV files from a directory, standardizes columns, adds IDs."""
#     print(f"\n--- Loading Metabolite Files from Directory: {metabolite_directory_path} ---")
#     processed_data = {} # Dictionary to store loaded DataFrames keyed by filename
#     if not os.path.isdir(metabolite_directory_path):
#         print(f"ERROR: Metabolite input directory not found: {metabolite_directory_path}")
#         return processed_data # Return empty dict
#     metabolite_file_dict = { f: os.path.join(metabolite_directory_path, f)
#                              for f in os.listdir(metabolite_directory_path) if f.lower().endswith('.csv') }
#     if not metabolite_file_dict:
#         print(f"Warning: No CSV files found in directory: {metabolite_directory_path}")
#         return processed_data # Return empty dict
#     print(f"Found {len(metabolite_file_dict)} potential CSV files to load.")

#     for filename, filepath in metabolite_file_dict.items():
#         print(f"\nProcessing: {filename}")
#         try:
#             df = pd.read_csv(filepath)
#             df.columns = [str(c).lower().strip() for c in df.columns] # Standardize cols
#             match = re.match(r'(MC\d+)', filename, re.IGNORECASE) or re.search(r'(MC\d+)', filename, re.IGNORECASE)
#             experiment_id = match.group(1).upper() if match else "UNKNOWN_EXPERIMENT"
#             if experiment_id == "UNKNOWN_EXPERIMENT": print(f"  Warning: No MC ID found in filename '{filename}'.")
#             df['experiment_id'] = experiment_id
#             df['source_filename'] = filename
#             processed_data[filename] = df
#         except Exception as e: print(f"  Error processing file {filename}: {e}")

#     if not processed_data: print("\nERROR: No metabolite files were successfully loaded.")
#     else: print(f"\nFinished loading metabolite files. {len(processed_data)} DataFrames loaded.")
#     return processed_data


# # --- Load RNA ---
# def load_rna_files(rna_directory_path, rows_are_genes=True):
#     """Loads multiple RNA seq count CSV files from a directory."""
#     print(f"\n--- Loading RNA Count Files from Directory: {rna_directory_path} ---")
#     loaded_data = {}
#     if not os.path.isdir(rna_directory_path):
#         print(f"ERROR: RNA input directory not found: {rna_directory_path}")
#         return loaded_data
#     rna_file_dict = { f: os.path.join(rna_directory_path, f)
#                       for f in os.listdir(rna_directory_path) if f.lower().endswith('.csv')}
#     if not rna_file_dict:
#         print(f"Warning: No CSV files found in directory: {rna_directory_path}")
#         return loaded_data
#     print(f"Found {len(rna_file_dict)} potential CSV files to load.")

#     for filename, filepath in rna_file_dict.items():
#         print(f"\nProcessing RNA file: {filename}")
#         try:
#             # Load as string initially, specify index column
#             # Use sep=None for auto-detection, engine='python' for flexibility
#             df_raw = pd.read_csv(filepath, index_col=0, sep=None, engine='python', dtype=str)
#             print(f"  Loaded '{filename}'. Initial shape: {df_raw.shape}")
#             loaded_data[filename] = df_raw
#         except Exception as e: print(f"  Error loading RNA count file '{filename}': {e}")

#     if not loaded_data: print("\nERROR: No RNA files were successfully loaded.")
#     else: print(f"\nFinished loading RNA files. {len(loaded_data)} DataFrames loaded.")
#     return loaded_data


# # --- Load RNA Metadata (Uses the local clean_index) ---
# def load_external_rna_metadata(rna_metadata_filepath):
#     """Loads external RNA metadata CSV, indexed by cleaned original_rna_sample_id."""
#     print(f"\n--- Loading External RNA Metadata: {os.path.basename(rna_metadata_filepath)} ---")

#     # 1. Check file existence
#     if not os.path.exists(rna_metadata_filepath):
#         print(f"ERROR: File not found at '{rna_metadata_filepath}'")
#         return None

#     # 2. Try loading, indexing, and cleaning
#     try:
#         # Assume first col is index ('original_rna_sample_id')
#         metadata_df = pd.read_csv(rna_metadata_filepath, index_col=0)

#         # Check if index was actually created
#         if not isinstance(metadata_df.index, pd.Index):
#              raise TypeError("Loaded data does not have a valid index after read_csv (check index_col=0).")

#         # *** Call the clean_index function defined WITHIN this file ***
#         metadata_df.index = clean_index(metadata_df.index)
#         metadata_df.index.name = 'original_rna_sample_id' # Ensure index name is set

#         print(f"  Successfully loaded and indexed. Shape: {metadata_df.shape}")
#         return metadata_df

#     except Exception as e: # Catch general exceptions during loading/indexing
#         print(f"ERROR loading or processing RNA metadata file '{os.path.basename(rna_metadata_filepath)}': {e}")
#         return None