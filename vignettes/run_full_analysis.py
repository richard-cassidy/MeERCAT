# Vignette
# Example Script using the MeERCAT package

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import meercat
from meercat import config
from meercat import load_data
from meercat import preprocess
from meercat import nmf
from meercat import visualize
from sklearn.preprocessing import LabelEncoder
import re


#from meercat import preprocess # Not using preprocessing in this example
INPUT_DATA_DIRECTORY = './input_data'

# --- Load Data ---
loaded_datasets = load_data.load_all_data(INPUT_DATA_DIRECTORY)

# --- Extract dataframes ---
rna_data = loaded_datasets.get('rna_data')
rna_metadata = loaded_datasets.get('rna_metadata')
experiment_metadata = loaded_datasets.get('experiment_metadata')
metabolite_data = loaded_datasets.get('metabolite_data')  

# --- Preprocess Data ---
transposed_rna_data = preprocess.transpose_rna_data(rna_data)
print("\n")
print("\n")
filtered_rna_data = preprocess.filter_rna_data(transposed_rna_data)
print("\n")
print("\n")
indexed_metabolite_data = preprocess.index_metabolite_data(metabolite_data)
print("\n")
print("\n")
combined_matched_data = preprocess.create_combined_dataframe(filtered_rna_data, indexed_metabolite_data )
print("\n")
print("\n")
print(combined_matched_data.columns)
print("\n")
print("\n")


#Focusing on metabolite data for now
metabs_scaled =nmf.scale_metabolites(combined_matched_data)

metabs_imputed = nmf.impute_for_nmf(metabs_scaled)

elbowplot= nmf.NMF_Elbow_Plot(metabs_imputed, range(2,41), nmf_params=None)

#This elbow plot showed that up to 17 components could be used. This is a high k value and the biology of the data that I know has shown that 7-10 types of metabolic processes are generlly changing in my system. With that in mind, I used k=10 for the analysis. 
print("\n")
print("\n")
print("\n")

nmf_params = {'max_iter': 500}
metab_W, metab_H = nmf.NMF_run(metabs_imputed, k=10, nmf_params=nmf_params)
print("\n")
print("\n")
print("\n")

#Able to run NMF and get H and W matrices 
test_name_H = metabs_scaled.columns.tolist()
test_name_W = metabs_scaled.index.tolist()
print(test_name_H)
print(test_name_W)


##Plotting/Data visualization and downstream analysis is still a work in progress. 
# Visualize with UMAP
visualize.visualize_nmf_with_umap(W= metab_W, H = metab_H, random_state=42,  show_labels = True, metric = "euclidean",title = "Test W values with name labels", feature_names = test_name_H, sample_names = test_name_W, matrix = "H") #set matrix value!
###Still working on coloring by sample groups/ k factors
