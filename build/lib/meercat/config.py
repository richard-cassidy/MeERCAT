# meercat_analyzer/meercat/config.py

# --- Default Paths ---
# Define default SUBDIRECTORY names. The main script will construct the full path.
# Use '.' to represent the current working directory or a user-defined base path.
DEFAULT_BASE_SAVE_PATH = './MeERCAT_Output/' # Default location for output relative to where script is run
DEFAULT_PROJECT_FOLDER = '' # Set to empty if output subfolders should be directly under base path,
                            # Or set a name like 'Analysis_Run' if you want BASE_PATH/Analysis_Run/Spearman...
DEFAULT_SPEARMAN_SUBFOLDER = 'Spearman_Analysis'
DEFAULT_NMF_SUBFOLDER = 'NMF_Analysis_Strategy2'
DEFAULT_PLOTS_SUBFOLDER = 'plots' # This will be relative within Spearman/NMF folders

# --- Default Filenames ---
METABOLITE_CLEANED_FILENAME = 'metabolite_data_cleaned_indexed.csv'
METADATA_EXTRACTED_FILENAME = 'metabolite_metadata_from_cleaning.csv'
RNA_COMBINED_FILENAME = 'rna_data_combined_raw.csv'
RNA_NORMALIZED_FILENAME = 'rna_data_normalized_log2.csv'
MATCHED_METABOLITE_FILENAME = 'metabolite_data_matched_aligned.csv'
MATCHED_RNA_FILENAME = 'rna_data_matched_aligned.csv'
CORRELATION_RAW_FILENAME = 'spearman_correlation_results_filtered.csv'
CORRELATION_ADJ_FILENAME = 'spearman_correlation_results_filtered_adj.csv'
NMF_H_FILENAME_TEMPLATE = 'NMF_Strategy2_H_matrix_k{k}.csv'
NMF_W_RNA_FILENAME_TEMPLATE = 'NMF_Strategy2_W_rna_matrix_k{k}.csv'
NMF_W_METAB_FILENAME_TEMPLATE = 'NMF_Strategy2_W_metab_matrix_k{k}.csv'

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
PLOT_HEATMAP_TOP_N_PAIRS = 200 # Increased from 100
PLOT_HEATMAP_MAX_FEATURES = 50 # Changed from 25
PLOT_N_TOP_POS_NEG_SCATTER = 10 # Changed from 15

# --- Metadata Columns (used for index creation, extraction, etc.) ---
# Adjust these based on your actual metadata/input files
METADATA_SAMPLE_COL = 'sample' # Column containing sample names before index setting
METADATA_ID_COLS = ['experiment_id', 'source_filename', 'group', 'ias_conc', 'rep', 'day', 'tracer']
METADATA_COMPOSITE_ID_COLS = ['experiment_id', 'ias_conc', 'day', 'rep'] # Used to build index
# Expected RNA metadata cols for separation before normalization
RNA_METADATA_COLS_EXPECTED = ['experiment_id', 'treatment_group', 'replicate', 'arsenic_concentration', 'days']

# --- Other ---
METABOLITE_ID_PREFIXES = ['ip_', 'hilica_', 'c30pos_']
METABOLITE_ID_SUFFIX_REGEX = r'(_ip|_hilica|_c30pos)(\.\d+)?$'
SAMPLE_FILTER_EXCLUDE_PATTERNS = ['qc', 'blank', 'kras']