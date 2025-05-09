# MeERCAT/src/meercat/preprocess.py
# process and combine RNA and metabolite to be ready for NMF
import pandas as pd
import pandas as pd
from typing import Dict






#Helper Function 1: Transpose RNA Data and indexes 
def transpose_rna_data(rna_data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Transposes the RNA DataFrames, making samples rows and genes columns.

    Args:
        rna_data: Dictionary of RNA DataFrames, keyed by experiment ID.
                  It expects each rna dataframe to have genes as rows, and sampleID's as its columns.
                  Where the column names become the index for the transposed DataFrames

    Returns:
        A dictionary of transposed RNA DataFrames.
    """
    transposed_rna_data = {}
    for exp_id, rna_df in rna_data.items():
        try:
            print(f"\n--- Transposing RNA Data for Experiment {exp_id} ---")
            transposed_df = rna_df.transpose()
            transposed_df = transposed_df.rename_axis('sample_ID')  # Make column to index

            transposed_rna_data[exp_id] = transposed_df
            print("\nTransposed output data")
            print(transposed_df)
        except Exception as e:
            print(f"  * Warning: Failed to transpose RNA data for experiment '{exp_id}': {e}")
            transposed_rna_data[exp_id] = rna_df # so it doesn't stop running
    return transposed_rna_data





#Helper Function 2: Filter RNA Data by count threshold
import pandas as pd

def filter_rna_at_threshold(rna_data, threshold=10, min_valid_fraction=0.9):
    """
    Filters RNA-seq data to remove rows that do not meet a minimum fraction of
    replicates with values above a threshold or not NA.

    Args:
        rna_data: Pandas DataFrame with RNA-seq counts.
        threshold: Minimum count value.
        min_valid_fraction: The rows must be at least over a percentage.

    Returns:
        Pandas DataFrame with filtered RNA-seq data.
    """
    # Number of Total Lines
    total_replicates = len(rna_data.columns)

    # Check with lines to keep
    valid_minimum_line = total_replicates*min_valid_fraction

    # Identify valid entries (not NA and above threshold)
    valid_values = (~rna_data.isnull()) & (rna_data >= threshold)

    # Count the number of valid values per row
    valid_count = valid_values.sum(axis=1)

    # Identify rows to keep: enough valid values
    keep_rows = valid_count >= valid_minimum_line

    # Apply the filter
    filtered_rna_data = rna_data[keep_rows]
    # print(f"\nFiltered RNA Data shape: {filtered_rna_data.shape}")
    # print(filtered_rna_data)
    return filtered_rna_data



def filter_rna_data(rna_data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Filters RNA DataFrames to keep only rows with index "ens_gene" or starting with 'MC00',
    sets the values from the row "ens_gene" as the new column names, and *removes* the
    "ens_gene" row.

    Args:
        rna_data: A dictionary of RNA DataFrames.

    Returns:
        A dictionary of filtered and renamed RNA DataFrames.
    """

    filtered_rna_data = {} # init
    for exp_id, rna_df in rna_data.items(): #for each
        try:
            print(f"\n--- Filtering RNA Data for Experiment {exp_id} ---")

            # Create a boolean mask to select rows starting with 'MC00'
            index_starts_with_mc00 = rna_df.index.str.startswith('MC00', na=False) # new matrix

            # Filter the DataFrame
            filtered_df = rna_df[index_starts_with_mc00] #

            # what is being renamed
            ens_gene_names = rna_df.loc['ens_gene'].values #
            # Rename the index
            filtered_df.columns = ens_gene_names #rename

            filtered_rna_data[exp_id] = filtered_df #assign as matrix, with values

            print(f"\nFiltered RNA Data for Experiment {exp_id}:")
            print(filtered_df)

        except Exception as e:
            print(f"  * Warning: Failed to filter RNA data for experiment '{exp_id}': {e}")
            filtered_rna_data[exp_id] = rna_df # and if it doesnt do that

    return filtered_rna_data





def index_metabolite_data(metabolite_data):
    indexed_metabolite_data = {}
    for exp_id, metabolite_df in metabolite_data.items():
        try:
            print(f"\nFiltering data in ex {exp_id} ")
            if 'sample_id' not in metabolite_df.columns:
                print(f"ERROR: {exp_id} .Returning what data I could")
                indexed_metabolite_data[exp_id] = metabolite_df
                continue

            metabolite_df = metabolite_df.set_index('sample_id')
            indexed_metabolite_data[exp_id] = metabolite_df

            print(f"Show  {exp_id}:")
            print(metabolite_df.head())
            #print(f"Show ID's {exp_id}: {metabolite_df.index.tolist()}") #print the


        except Exception as e:
            print(f"ERROR DATA IN {exp_id}: {e}")
            indexed_metabolite_data[exp_id] = metabolite_df

    return indexed_metabolite_data





#Helper Function 3: Drop columns with high NaN

def drop_cols_with_high_na(df: pd.DataFrame, na_threshold: float = 0.5) -> pd.DataFrame:
    """
    Removes columns from a Pandas DataFrame if the percentage of NaN values
    in that column is greater than a specified threshold.

    Args:
        df (pd.DataFrame): The input DataFrame.
        na_threshold (float): The threshold for the percentage of NaN values.
                              Columns with a higher percentage of NaNs will be dropped.
                              Defaults to 0.9 (90%).

    Returns:
        pd.DataFrame: A new DataFrame with columns exceeding the NaN threshold removed.
    """

    # Calculate the percentage of NaN values in each column
    print(len(df))
    print(df)
    na_percentages = df.isnull().sum() / len(df)
    #print(na_percentages)
    # Identify columns exceeding the NaN threshold
    cols_to_drop = na_percentages[na_percentages > na_threshold].index
    # len(cols_to_drop)
    # print(len(cols_to_drop))
    # print(f"Columns to drop: {cols_to_drop}")
    # Drop columns
    df = df.drop(cols_to_drop, axis=1)

    return df









def create_combined_dataframe(rna_data: Dict[str, pd.DataFrame], metabolite_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Creates a single DataFrame containing both RNA and metabolite data for matching indices.
    Removes non-numeric columns, then removes columns with high missing values (NaNs),
    and finally imputes the remaining missing values with the column median *within each experiment*.

    Args:
        rna_data: Dictionary of RNA DataFrames (keyed by experiment ID). Assumes data has index.
        metabolite_data: Dictionary of Metabolite DataFrames (keyed by experiment ID). Assumes data has index.

    Returns:
        A single Pandas DataFrame with combined RNA and metabolite data for matching indices.
        Returns an empty DataFrame if rna_data or metabolite_data is empty,
        or if no shared indices are found.
    """

    combined_data = []  # to generate combined matrices

    try:
        exp_ids_in_both = set(rna_data.keys()) & set(metabolite_data.keys())  # Find experiments with both

        if not exp_ids_in_both:  # if non are in this
            print("Warning: No shared experiments")
            return pd.DataFrame()  # to quit and return empty

        for exp_id in exp_ids_in_both:
            rna_df = rna_data[exp_id]
            metabolite_df = metabolite_data[exp_id]  # read in now and only to

            # Get shared index
            shared_index = rna_df.index.intersection(metabolite_df.index)
            if shared_index.empty:
                print(f"No shared indices between RNA and metabolite data for experiment {exp_id}. Skipping")
                continue

            # Combine using a left join on the shared index, keeping RNA data rows
            tempCombined = pd.merge(rna_df.loc[shared_index], metabolite_df.loc[shared_index], left_index=True, right_index=True, how='left')
            tempCombined["experiment_id"] = exp_id
            print("\n--- tempCombined BEFORE filtering ---")
            print(tempCombined.head()) # show the first few rows
            print(tempCombined.isnull().sum()) # check the number of NaNs per column
            print(tempCombined.dtypes) # check datatypes.

            # ***ADDED: Remove Non-Numeric Columns***
            numeric_cols = tempCombined.select_dtypes(include=['number']).columns #select columns that include only numbers - int/float/etc
            tempCombined = tempCombined[numeric_cols] #select those columns

            # Drop columns with high NaN values *before* imputation
            tempCombined = drop_cols_with_high_na(tempCombined)
            print("\n--- tempCombined AFTER filtering ---")
            print(tempCombined.head()) # show the first few rows
            print(tempCombined.isnull().sum()) # check the number of NaNs per column
            print(tempCombined.dtypes) # check datatypes.           
            # Impute missing values with median, calculated *within* each column and DataFrame
            tempCombined = tempCombined.fillna(tempCombined.median())

            combined_data.append(tempCombined)  # Collect combined experiments

        # Concatenate all experiment DataFrames
        final_df = pd.concat(combined_data, axis=0)

        # Drop columns with high NaN across *all* combined data (if any remain after per-experiment processing)
        final_df = drop_cols_with_high_na(final_df)

        print(f"\nCombined Data shape {final_df.shape} has compiled - it's there so now you see if its right. Best of luck : ")
        return final_df

    except Exception as e:
        print(f"There was a great mistake in making it one: {e}")
        return pd.DataFrame()