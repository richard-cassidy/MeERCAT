# MeERCAT/src/meercat/load_data.py
# Load Data Function


import pandas as pd
import os
from typing import Dict, Optional
from meercat import config
import pandas as pd
import os
from typing import Dict, Optional
import pandas as pd
import anndata as ad
import numpy as np


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
                        try: 
                          rna_data[exp_id] = pd.read_csv(file_path)
                          if rna_data[exp_id] is not None:
                            rna_data[exp_id].columns = rna_data[exp_id].columns.str.lower()
                            # INSPECTION POINT 1: Print RNA data shape and columns immediately after loading
                            print(f"\n--- RNA Data for Experiment {exp_id} (After Loading) ---")
                            print(f"Shape: {rna_data[exp_id].shape}")
                            print(f"Columns: {rna_data[exp_id].columns.tolist()}")
                        except Exception as e: print(f"  * Warning: Failed to load RNA file {filename}: {e}")
                elif filename.endswith('_metabolites.csv'):
                    exp_id = filename.removesuffix('_metabolites.csv')
                    if exp_id:
                       try: 
                          metabolite_data[exp_id] = pd.read_csv(file_path)
                          if metabolite_data[exp_id] is not None:
                            metabolite_data[exp_id].columns = metabolite_data[exp_id].columns.str.lower()
                       except Exception as e: print(f"  * Warning: Failed to load Metabolite file {filename}: {e}")
        print(f"   loaded RNA and metabolite files")
        print(f"----------------------------\n\n")

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

        # INSPECTION POINT 2: Print RNA metadata shape, columns, and unique values
        print("\n--- RNA Metadata (After Loading) ---")
        print(f"Shape: {metadata['rna'].shape}")
        print(f"Columns: {metadata['rna'].columns.tolist()}")
        print(f"Unique original_rna_sample_id: {metadata['rna']['original_rna_sample_id'].unique().tolist()}")
        print(f"Unique sample_id: {metadata['rna']['sample_ID'].unique().tolist()}")

    except FileNotFoundError:
        print(f"  * Warning: RNA metadata file not found: {rna_meta_path}")
    except Exception as e:
        print(f"  * Warning: Failed to load RNA metadata file {rna_metadata_filename}: {e}")
    return metadata

import pandas as pd


def rename_rna_data_columns(
    rna_data: dict[str, pd.DataFrame], rna_metadata: pd.DataFrame
) -> dict[str, pd.DataFrame]:
    """Renames columns in RNA data DataFrames using the RNA metadata.
    Handles case sensitivity, whitespace, and missing values.
    """

    if rna_metadata is None or rna_metadata.empty:
        print(
            "RNA metadata is missing or empty. Returning RNA data with original"
            " column names."
        )
        return rna_data

    renamed_rna_data = {}
    for exp_id, rna_df in rna_data.items():
        try:
            # INSPECTION POINT 3: Print RNA data columns BEFORE renaming in this loop
            print(f"\n--- RNA Data for Experiment {exp_id} (Before Renaming) ---")
            print(f"Columns: {rna_df.columns.tolist()}")

            # Convert relevant columns to string type and strip whitespace:
            rna_df.columns = rna_df.columns.astype(str).str.strip()
            rna_metadata['original_rna_sample_id'] = rna_metadata['original_rna_sample_id'].astype(str).str.strip()
            rna_metadata['sample_ID'] = rna_metadata['sample_ID'].astype(str).str.strip()

            # Convert to lowercase for case-insensitive matching
            rna_metadata['original_rna_sample_id'] = rna_metadata['original_rna_sample_id'].str.lower()
            rna_df.columns = rna_df.columns.str.lower()

            # Filter RNA metadata, remove duplicates, and handle NaN values
            filtered_rna_metadata = rna_metadata[
                ["original_rna_sample_id", "sample_ID"]
            ].drop_duplicates().dropna(subset=['original_rna_sample_id', 'sample_ID'])  # Drop rows with NaN in relevant columns

            # INSPECTION POINT 4: Print filtered metadata before mapping
            print("\n--- Filtered Metadata ---")
            print(filtered_rna_metadata.head())

            # Create mapping dictionary
            id_mapping = dict(
                zip(
                    filtered_rna_metadata["original_rna_sample_id"],
                    filtered_rna_metadata["sample_ID"],
                )
            )

            # INSPECTION POINT 5: Print the ID mapping dictionary
            print("\n--- ID Mapping Dictionary ---")
            print(id_mapping)

            # Rename columns in the RNA DataFrame
            new_columns = []
            for col in rna_df.columns:
                if col in id_mapping:
                    new_columns.append(id_mapping[col])
                else:
                    new_columns.append(col)  # Keep original if not in metadata

            renamed_rna_data[exp_id] = rna_df.copy()
            renamed_rna_data[exp_id].columns = new_columns

        except Exception as e:
            print(
                f"  * Warning: Failed to rename RNA data columns for experiment"
                f" '{exp_id}': {e}.  Returning original DataFrame."
            )
            renamed_rna_data[exp_id] = rna_df
            print(f"Experiment {exp_id}: id_mapping = {id_mapping}")
        
    return renamed_rna_data


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

    # --- Rename RNA data columns ---
    all_data["rna_data"] = rename_rna_data_columns(all_data["rna_data"], all_data["rna_metadata"])


    # --- Final Summary ---
    print("\nSuccessfully Loaded")
    if successfully_loaded_list:
        for item in successfully_loaded_list:
            print(f"  - {item}")
    else:
        print("  * Warning: No data components were successfully loaded.")
    print("----------------------------\n\n")


    return all_data