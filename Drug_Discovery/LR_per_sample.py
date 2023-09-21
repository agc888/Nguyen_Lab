#!/usr/bin/env python
# coding: utf-8

# In[1]:


import stlearn as st
from pathlib import Path
st.settings.set_figure_params(dpi=180)
import pandas as pd
import math
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import pickle
import sys



# # Import Data

with open('/QRISdata/Q1851/Andrew_C/Breast_Cancer_Drug/cancer_only.pickle', 'rb') as f:
    data_dict = pickle.load(f)


samples = list(data_dict.keys())
i = int(sys.argv[1])
sample = samples[i]



DATA_PATH = "/QRISdata/Q2051/Pfizer/Visium/RAW_DATA/Pfizer/"
OUT_PATH = "/QRISdata/Q1851/Andrew_C/Breast_Cancer_Drug/"

data = data_dict[sample]

## Make file save directory

lr_file_path = OUT_PATH + sample+"/"
stlearn_dir = Path(lr_file_path)

# Check if the directory already exists
if not stlearn_dir.exists():
    # Create the directory
    stlearn_dir.mkdir()


### RUN STLEARN L-R Analysis ###

lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')

st.tl.cci.run(data, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=0, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=None, # Number of CPUs for parallel. If None, detects & use all available.
                  ) 


st.tl.cci.adj_pvals(data, correct_axis='spot',pval_adj_cutoff=0.05, adj_method='fdr_bh')


# Save Output

data.uns["lrfeatures"] = data.uns["lrfeatures"].astype(str)
data.write_h5ad(lr_file_path+sample+"_lr_data.h5ad")

