
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import NMF
from sklearn.impute import KNNImputer

import matplotlib.pyplot as plt












def scale_rna(rna_data):
    """Scales RNA-seq data (CPM normalization, log1p transformation, and standard scaling).
    Selects only numerical columns with names starting with "ENSG".
    Removes rows that are not numeric.
    """
    # Select columns with names starting with "ENSG"
    print("\n\n RNA data")
    print(rna_data)
    print("\n\n")
    print("\n\n RNA data columns")
    print(rna_data.columns)
    print("\n\n")
    ens_cols = [col for col in rna_data.columns if str(col).startswith("ENSG")]
    rna_data = rna_data[ens_cols]

    # Identify rows with any non-numeric values
    numeric_check = rna_data.applymap(lambda x: isinstance(x, (int, float)))
    non_numeric_rows = numeric_check[~numeric_check.all(axis=1)].index

    # Drop the non-numeric rows
    rna_data = rna_data.drop(non_numeric_rows)

    print("\n\nRNA Data:")
    print(rna_data)

    # CPM normalization, log1p transformation, and standard scaling
    cpm_data = rna_data.apply(lambda x: x / x.sum() * 1e6, axis=0)
    log_transformed = np.log1p(cpm_data)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(log_transformed)
    return scaled_data














def scale_metabolites(metabolite_data: pd.DataFrame):
    """Scales metabolite data (mean normalization, log1p transformation, and standard scaling),
    handling NaNs correctly and ensuring non-negativity.
    Selects only numerical columns, excluding those that start with "ENS".
    """
    # Select numerical columns, excluding those that start with "ENS"
    metabolite_cols = [col for col in metabolite_data.columns if not str(col).startswith("ENS") and pd.api.types.is_numeric_dtype(metabolite_data[col])]
    metabolite_data = metabolite_data[metabolite_cols].copy() # Create a copy
    # plt.hist(metabolite_data.values.flatten(), bins=50, color='skyblue', edgecolor='black')
    # print("\n\nMetabolite Data:")
    # print(metabolite_data.head())  # Print a preview

    # Mean normalization, handling NaNs: mean will ignore NaNs by default
    mean_normalized = metabolite_data.apply(lambda x: x / x.mean(), axis=0)
    #Log Transform
    #print((mean_normalized < 0).any().any()) # will print "True" if there are any such values


    #Multiply by ~1000 or so 
    mean_norm_x10000= mean_normalized * 10000
   
    # print((mean_norm_x10000 < 0).any().any()) # will print "True" if there are any such values

    # Log1p transformation, handling NaNs: log1p handles NaNs gracefully
    log_transformed = np.log1p(mean_norm_x10000)
    log_transformed = log_transformed-6.5
    log_transformed[log_transformed < 0] = 0
    plt.hist(log_transformed.values.flatten(), bins=50, color='skyblue', edgecolor='black')
    plt.title('Distribution of Scaled Metabolite Data')
    plt.xlabel('Scaled Value')
    plt.ylabel('Frequency')

    #save plot as a PNG file
    plt.savefig('scaled_metabolite_distribution.png')
    #plt.show()

    return log_transformed










def impute_for_nmf(data, threshold=0.9, n_neighbors=5):
    """
    Prepares data for NMF. Drops columns with high missingness, then imputes
    remaining missing values (NaNs) using KNN imputation.

    Args:
        data: The input DataFrame.
        threshold: The threshold for the percentage of NaN values. Columns with a
                   higher percentage of NaNs will be dropped. Defaults to 0.9 (90%).
        n_neighbors: The number of neighbors to use for KNN imputation (default: 5).

    Returns:
        matrix: The data matrix, with high-NaN columns dropped and remaining NaNs
                imputed using KNN, as a NumPy array.
    """

    # Identify numeric columns
    numeric_cols = data.select_dtypes(include=np.number).columns

    # Keep only numeric columns
    data = data[numeric_cols]

    # Identify columns with too many missing values
    missing = data.isnull().sum() / len(data)
    drop = missing[missing > threshold].index
    data = data.drop(columns=drop)

    # Convert to float type if not already (important for KNNImputer)
    data = data.astype(float)

    # Initialize KNNImputer
    imputer = KNNImputer(n_neighbors=n_neighbors)

    # Impute missing values
    matrix = imputer.fit_transform(data)

    # Set values less than 0 to 0
    matrix[matrix < 0] = 0

    return matrix




















def NMF_Elbow_Plot(data: np.ndarray, k_range: range, nmf_params: dict = None):
    """
    Generates an Elbow plot to aid in selecting the optimal number of components (k)
    for Non-negative Matrix Factorization (NMF). The plot visualizes the reconstruction
    error for different values of k.

    Args:
        data: The input data as a NumPy array (output from prep_imputed_nmf).
        k_range: A range object specifying the range of k values to test (e.g., range(2, 11)).
        nmf_params: A dictionary with NMF hyperparameter options (e.g., n_init, max_iter).
    """

    wcss = []  # Within-cluster sum of squares (inertia)
    nmf_params = nmf_params or {}  # Use default empty dict if nmf_params is None
    print(f"Input data shape: {data.shape}")
    print(f"Data contains NaN values: {np.isnan(data).any()}")
    print(f"Data contains infinite values: {np.isinf(data).any()}")
    print(f"Minimum data value: {np.min(data)}")
    print(f"Maximum data value: {np.max(data)}")
    print(f"Average data value: {np.average(data)}")
    if data.size > 0:  # Avoid error if data is empty
        print(f"Data contains negative values: {(data < 0).any()}")

        
    for k in k_range:
        print(f"Fitting NMF with k={k}")
        print(f"Shape: {data.shape}")
        # Initialize NMF model with specified hyperparameters
        nmf_model = NMF(n_components=k, random_state=42, **nmf_params)

        # Fit NMF model to the data
        nmf_model.fit(data)
        print(f"K={k} - NMF fit successful, err {nmf_model.reconstruction_err_}")
        # Calculate reconstruction error (inertia or WCSS)
        reconstruction_error = nmf_model.reconstruction_err_
        wcss.append(reconstruction_error)
        



    # Plot the Elbow curve
    plt.figure(figsize=(10, 6))
    plt.plot(k_range, wcss, marker='o', linestyle='-', color='skyblue', linewidth=2, markersize=8, markerfacecolor='royalblue')
    plt.title('Elbow Method for Optimal k')
    plt.xlabel('Number of Components (k)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS) / Reconstruction Error')
    plt.grid(True)  # Add a grid for better readability
    plt.xticks(k_range)  # Ensure all k values are visible on the x-axis
    plt.savefig("elbow_curve.png")  # saves curve as a png
    #plt.show()  # view plot

    print("Elbow Plot saved as elbow_curve.png. Please visually inspect the plot to determine the optimal k.")



















def NMF_run(data: np.ndarray, k: int, nmf_params: dict = None):
    """
    Runs Non-negative Matrix Factorization (NMF) on the input data with specified
    number of components (k) and hyperparameters.

    Args:
        data: The input data as a NumPy array (output from prep_imputed_nmf).
        k: The number of components for NMF.
        nmf_params: A dictionary with NMF hyperparameter options (e.g., n_init, max_iter).

    Returns:
        W: The basis matrix (features).
        H: The coefficient matrix (weights).
    """
    nmf_params = nmf_params or {}  # Use default empty dict if nmf_params is None
    print(f"Input data shape: {data.shape}")
    print(f"Data contains NaN values: {np.isnan(data).any()}")
    print(f"Data contains infinite values: {np.isinf(data).any()}")
    print(f"Minimum data value: {np.min(data)}")
    print(f"Maximum data value: {np.max(data)}")
    print(f"Average data value: {np.average(data)}")
    if data.size > 0:  # Avoid error if data is empty
        print(f"Data contains negative values: {(data < 0).any()}")

        
    # Initialize NMF model with specified hyperparameters
    nmf_model = NMF(n_components=k, random_state=42, **nmf_params)

    # Fit NMF model to the data
    W = nmf_model.fit_transform(data)
    print(W)
    print(f"Shape of W: {W.shape}")

    H = nmf_model.components_
    print(H)
    print(f"Shape of H: {H.shape}")
    
    return W, H





