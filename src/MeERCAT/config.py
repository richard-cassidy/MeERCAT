# MeERCAT/src/meercat/config.py



# --- Default Paths ---
# Define default SUBDIRECTORY names. The main script will construct the full path.
# Use '.' to represent the current working directory or a user-defined base path.
DEFAULT_BASE_SAVE_PATH = './meercat_output/' # Default location for output relative to where script is run
DEFAULT_PROJECT_FOLDER = '' # Set to empty if output subfolders should be directly under base path,
DEFAULT_SPEARMAN_SUBFOLDER = 'correlation_analysis'
DEFAULT_NMF_SUBFOLDER = 'NMF_Analysis'
DEFAULT_PLOTS_SUBFOLDER = 'plots' # This will be relative within Spearman/NMF folders

# --- Default Filenames ---
# Input Data Filenames
##### Assumed relative to input_data folder set in main script.

RNA_METADATA_FILENAME = 'rna_metadata.csv' # Metadata for RNA samples
EXPERIMENT_MATADATA_FILENAME = 'experiment_metadata.csv' # General metadata for all experiments
RNA_EXTERNAL_METADATA_COLS = ['sample_ID', 'ias_conc', 'day', 'rep']


# Output Data Filenames
##### Assumed relative to base output path

#Processed metabolite data
METABOLITE_CLEANED_FILENAME = 'metabolite_data_cleaned_indexed.csv'
METABOLITE_METADATA_EXTRACTED_FILENAME = 'metabolite_metadata_from_cleaning.csv'

#Processed RNA data
RNA_CLEANED_AND_NORMALIZED_FILENAME = 'rna_data_cleaned_and_normalized.csv'

#Merged RNA and Metabolite data from cleaned datasets
RNA_METABOLITES_ALL_MERGED_FILENAME = 'rna_metabolite_all_merged.csv'
RNA_METABOLITES_ONLY_MATCHED_FILENAME = 'rna_metabolite_matched.csv'

#Correlation files
CORRELATION_RAW_FILENAME = 'spearman_correlation_results_filtered.csv'
CORRELATION_ADJ_FILENAME = 'spearman_correlation_results_filtered_adj.csv'

#NMF files
NMF_H_FILENAME_TEMPLATE = 'NMF_H_matrix_k{k}.csv'
NMF_W_RNA_FILENAME_TEMPLATE = 'NMF_W_rna_matrix_k{k}.csv'
NMF_W_METAB_FILENAME_TEMPLATE = 'NMF_W_metab_matrix_k{k}.csv'

# --- Filtering Parameters ---
VAR_FILTER_TOP_N_GENES = 500
VAR_FILTER_TOP_N_METABOLITES = 150
RNA_FILTER_MIN_COUNTS = 10
RNA_FILTER_MIN_SAMPLES_FRAC = 0.2

# --- NMF Parameters ---
NMF_DEFAULT_K = 5
NMF_DEFAULT_MAX_ITER = 5000
NMF_DEFAULT_RANDOM_STATE = 42

# --- Plotting Parameters ---
PLOT_PADJ_THRESHOLD = 0.05
PLOT_RHO_THRESHOLD = 0.6
PLOT_N_TOP_VOLCANO_LABELS = 10
PLOT_HEATMAP_TOP_N_PAIRS = 200
PLOT_HEATMAP_MAX_FEATURES = 50
PLOT_N_TOP_POS_NEG_SCATTER = 10

# --- Metadata Columns (used for index creation, extraction, etc.) ---
# Adjust these based on your actual metadata/input files
METADATA_SAMPLE_COL = 'sample' # Column containing sample names before index setting (in raw metabolite data)
# List of potential metadata columns that might be extracted or present
METADATA_ID_COLS = ['experiment_id', 'source_filename', 'group', 'ias_conc', 'rep', 'day', 'tracer']
# Columns REQUIRED to build the standard composite ID
METADATA_COMPOSITE_ID_COLS = ['experiment_id', 'ias_conc', 'day', 'rep']
# Columns expected in the EXTERNAL RNA metadata file (besides original_rna_sample_id index)
# These should align with METADATA_COMPOSITE_ID_COLS for consistent ID creation
#RNA_EXTERNAL_METADATA_COLS = ['experiment_id', 'ias_conc', 'day', 'rep'] # Make sure these match Composite ID cols
# Columns expected in the combined RNA dataframe BEFORE normalization (for separation)
RNA_METADATA_COLS_EXPECTED = ['experiment_id', 'treatment_group', 'replicate', 'arsenic_concentration', 'days'] # Used in normalize function - needs review based on actual processing

# --- Other ---
METABOLITE_ID_PREFIXES = ['ip_', 'hilica_', 'c30pos_']
METABOLITE_ID_SUFFIX_REGEX = r'(_ip|_hilica|_c30pos)(\.\d+)?$'
SAMPLE_FILTER_EXCLUDE_PATTERNS = ['qc', 'blank', 'kras']