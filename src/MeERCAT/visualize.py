# mercat_analyzer/mercat_analyzer/visualize.py

import umap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re




def visualize_nmf_with_umap(
    W: np.ndarray, #num_samples x components
    H: np.ndarray, #num_features x components
    sample_names: list = None,
    feature_names: list = None,
    sample_group_assignments: list = None,
    title: str = "UMAP Visualization of NMF Results",
    show_labels = False,
    matrix: str = "H", # W matrix or H matrix
    n_components: int = 2, # UMAP dimensions
    random_state: int = 42,
    n_neighbors:int = 10,
    min_dist:float = 0.1,
    metric:str = 'euclidean'
) -> None:
    """
    Visualizes NMF results (W or H matrix) using UMAP.

    Args:
        W: The basis matrix (samples x components).
        H: The coefficient matrix (components x features).
        sample_names: List of sample names for the rows of W (optional).
        feature_names: List of feature names for the columns of H (optional).
        sample_group_assignments: List of group assignments for samples (optional).
            If provided, samples will be colored by group.
        title: Plot title.
        show_labels: display sample/feature name on scatterplot, or sample/feature group assignments

    Raises:
        ValueError: If sample_group_assignments is provided but its length does not
            match the number of samples (rows in W).
        ValueError: If both sample_group_assignments and feature_names are provided
            as they are mutually exclusive data.
    """
    #Validate Data

    if matrix != 'W' and matrix != "H": #check is good data
        raise ValueError(f"If grouping/colors are selected, must be string: 'W' or 'H' NOT {matrix}")
    if sample_group_assignments is not None and len(sample_group_assignments) != W.shape[0]:
        raise ValueError(
            "Length of sample_group_assignments does not match the number of samples"
        )
    if sample_group_assignments is not None and feature_names is not None:
        raise ValueError(
            "can display sample groups or feature names, these are mutually exclusive parameters - please pick one"
        )
        
    if matrix == "H":
        embedding = umap.UMAP(n_components=n_components, random_state=random_state, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric).fit_transform(H.T) # H values are passed on
    else:
        embedding = umap.UMAP(n_components=n_components, random_state=random_state, n_neighbors = n_neighbors, min_dist = min_dist, metric = metric).fit_transform(W) # H values are passed on
  

    # Generate base visualization of the UMAP

    if show_labels == True:
        if matrix == "H": # for matrix that's H
             colors_assigned = list(feature_names)
        else: # otherwise its matrix W
             colors_assigned = list(sample_names)
        fig, ax = plt.subplots()

        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.5)
        plt.xlabel("UMAP Dimension 1")
        plt.ylabel("UMAP Dimension 2")
        plt.title(title)
        if colors_assigned != None:

           for i in range(embedding.shape[0]):
            ax.text(embedding[i, 0], embedding[i, 1], str(colors_assigned[i]), fontsize = 3)
        else:
           print ("Check names of Features to show are set correctly - or no label chosen")

    elif sample_group_assignments is not None: #otherwise use the assigned data to sample data

       fig, ax = plt.subplots()

       plt.scatter(embedding[:, 0], embedding[:, 1], c=sample_group_assignments, alpha=0.5)
       plt.xlabel("UMAP Dimension 1")
       plt.ylabel("UMAP Dimension 2")
       plt.title(title)

    else:
       fig, ax = plt.subplots()

       plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.5)
       plt.xlabel("UMAP Dimension 1")
       plt.ylabel("UMAP Dimension 2")
       plt.title(title)
   # Show or return the UMAP:
    plt.savefig("umap_plot.png", dpi=300)  # Save the plot as a PNG file













from sklearn.preprocessing import LabelEncoder
import re

def extract_arsenic_concentration(sample_name):
  """Extracts the arsenic concentration from the sample name."""
  match = re.search(r'_(\d+)_', sample_name)  # Find the number between underscores
  if match:
    return match.group(1)  # Return the captured group (the concentration)
  else:
    return None  # Or handle the case where the pattern isn't found

