;; This buffer is for text that is not saved, and for Lisp evaluation.
;; To create a file, visit it with <open> and enter text in its buffer.

import pandas as pd
import scimap as sm
import anndata as ad

def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.items():
        print("    %s: %s" % (key, val))

def print_all_keys(h5_file):
    def print_key(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(name)

    h5_file.visititems(print_key)

f = h5py.File('output.h5','r');    
with h5py.File('example.h5ad', 'r') as f:
    #print_all_keys(f)
    f.visititems(print_attrs)
    X = f['X'][:]
    rvar = f['var'][:]
    # Decode byte strings to regular strings
    rvar = [var.decode('utf-8') for var in rvar]

# Convert obs and var to pandas DataFrames
obs = pd.DataFrame(index=robs)
var = pd.DataFrame(index=rvar)

obs = pd.DataFrame(robs)
obs = obs.applymap(lambda x: x.decode('utf-8'))
obs.index = range(obs.shape[0])  # this sets the index to be from 0 to n-1
var = pd.DataFrame(rvar)
var = var.applymap(lambda x: x.decode('utf-8'))
var.index = range(var.shape[0])  # this sets the index to be from 0 to n-1

# Construct AnnData object
adata = ad.read_h5ad('output.h5');
    
with h5py.File('output.h5', 'r') as f:
    # Now you can access data in the file
    for key in f.keys():
        print(key)  # Prints the names of all top-level groups

aj = ad.read_h5ad('example.h5ad')        
