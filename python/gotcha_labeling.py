# -*- coding: utf-8 -*-
"""
@author: Sanjay Kottapalli (svk4001@med.cornell.edu)
@date: 08/11/2022
Landau Lab, WCM

This script executes functions for the genotyping of cells, 
starting from wildtype and MUT read counts.
"""

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import timeit
from sklearn.neighbors import KernelDensity
from collections import Counter
from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn.neighbors import KNeighborsClassifier
from scipy.stats import gmean

# Set styles for plots
plt.rcParams['figure.dpi'] = 600
sns.set(font_scale=1.2)
sns.set_style("ticks")
tab20 = plt.get_cmap('tab20')

def GotchaLabeling(df=None, merged_genotype_file="", saturation=False):
    time1 = timeit.default_timer()
    print("Reading in file.")
    typing = read_data(df, merged_genotype_file)
    
    if typing is None:
        return None
    
    parent_dir = os.path.dirname(merged_genotype_file) if merged_genotype_file else os.getcwd()
    sample_dir = os.path.join(parent_dir, '05outs/')
    os.makedirs(sample_dir, exist_ok=True)
    
    for i in ["wildtype", "mutant"]:
        print(f"Noise correcting {i} read counts.")
        typing, min_val, thresh, _ = noise_correct(typing, i, sample_dir)
        
        if typing is None:
            return None
    
    typing = quadrant_genotype(typing, min_val, min_val)
    typing = KNN_cluster(typing, min_val, min_val, 5, sample_dir)
    
    typing.to_csv(os.path.join(sample_dir, '_genotype_labels.csv'))
    
    print("All analysis complete!")
    time2 = timeit.default_timer()
    print(f"Total time to execute: {time2 - time1}")
    
    return typing

def read_data(df=None, merged_genotype_file="", barcode_column="barcode"):
    if df is None:
        cell_line = pd.read_csv(merged_genotype_file, index_col=0, sep=",")
    else:
        cell_line = df
        
    cell_line.fillna(0.0, inplace=True)

    print("Data loaded:")
    print(cell_line.head())

    try:
        genotyping = pd.DataFrame(index=cell_line.index)
        genotyping['wildtype_count'] = cell_line['wildtype_count']
        genotyping['mutant_count'] = cell_line['mutant_count']
    except KeyError as e:
        print(f"Missing column: {e}")
        return None

    index = genotyping[['wildtype_count', 'mutant_count']].dropna().index
    genotyping = genotyping.loc[index, :]
    
    print("Number of cells: " + str(genotyping.shape[0]))
    
    return genotyping

def noise_correct(typing, feature="", sample_dir=""):
    np.random.seed(0)
    pseudocount = 1
    X = typing[feature + '_count'] + pseudocount
    logged_counts = np.log(X)

    typing['transf_{}'.format(feature)] = logged_counts

    plt.hist(logged_counts, density=True, bins=50)
    plt.title(f"{feature} counts")
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.show()
    plt.clf()

    bw = 0.1
    kde = KernelDensity(bandwidth=bw)
    kde.fit(typing['transf_{}'.format(feature)].values.reshape(-1, 1))

    x_bin = np.histogram(typing['transf_{}'.format(feature)], bins=50)[1]
    kde_x = np.linspace(min(x_bin) - 0.5, max(x_bin) + 0.5, 10000)
    kde_smooth = np.exp(kde.score_samples(kde_x.reshape(-1, 1)))

    plt.hist(typing['transf_{}'.format(feature)], density=True, bins=50)
    plt.plot(kde_x, kde_smooth, color='red')
    plt.title('Initial KDE')
    plt.ylabel("Probability")
    plt.xlabel("Log(counts+1)")
    plt.savefig(os.path.join(sample_dir, "{}_kde_initial.pdf".format(feature)), dpi=500, bbox_inches="tight")
    plt.show()
    plt.clf()

    noise_values = logged_counts[logged_counts < min(kde_x)]
    noise_values = noise_values.dropna()

    if len(noise_values) == 0:
        print("No noise values found. Exiting function.")
        return typing, None, None, None

    kde_noise = KernelDensity(bandwidth=bw)
    kde_noise.fit(noise_values.values.reshape(-1, 1))
    noise_smooth = np.exp(kde_noise.score_samples(kde_x.reshape(-1, 1)))

    signal_values = logged_counts[logged_counts >= min(kde_x)]
    signal_values = signal_values.dropna()

    if len(signal_values) == 0:
        print("No signal values found. Exiting function.")
        return typing, None, None, None

    kde_signal = KernelDensity(bandwidth=bw)
    kde_signal.fit(signal_values.values.reshape(-1, 1))
    signal_smooth = np.exp(kde_signal.score_samples(kde_x.reshape(-1, 1)))

    # Additional analysis...
    
    return typing, None, None, None  # Adjust return values as needed

def quadrant_genotype(typing, wt_min, mut_min):
    typing['quadrant_class'] = 'NA'
    for i in typing.index:
        if typing.loc[i, 'transf_wildtype'] < wt_min:
            if typing.loc[i, 'transf_mutant'] < mut_min:
                typing.loc[i, 'quadrant_class'] = 'NA'
            else:
                typing.loc[i, 'quadrant_class'] = 'MUT'
        else:
            if typing.loc[i, 'transf_mutant'] < mut_min:
                typing.loc[i, 'quadrant_class'] = 'WT'
            else:
                typing.loc[i, 'quadrant_class'] = 'HET'
    
    return typing

def KNN_cluster(typing, wt_min, mut_min, knn_window=0.05, sample_dir=""):
    typing['clusters'] = typing['quadrant_class']
    data = typing[['transf_wildtype', 'transf_mutant']].values
    range_wt = max(data[:, 0]) - min(data[:, 0])
    range_mut = max(data[:, 1]) - min(data[:, 1])
    
    indices1 = set(np.where(data[:, 0] > wt_min - knn_window * range_wt)[0])
    indices1 = indices1.intersection(set(np.where(data[:, 0] <= wt_min + knn_window * range_wt)[0]))
    indices2 = set(np.where(data[:, 1] > mut_min - knn_window * range_mut)[0])
    indices2 = indices2.intersection(set(np.where(data[:, 1] <= mut_min + knn_window * range_mut)[0]))
    indices = indices1.union(indices2)
    
    indices = np.array(list(indices))
    indices = typing.index[indices]
    
    typing.loc[indices, 'clusters'] = -1
    
    n = list(dict(Counter(typing['quadrant_class'].values)).values())
    print(n)
    n_neighbors = round(np.sqrt(gmean(n)))
    print("Nearest neighbours: " + str(n_neighbors))
    
    print("Remove uncertain labels and re-label with self-training.")
    STC = SelfTrainingClassifier(KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance'),
                                  criterion='k_best', k_best=1, max_iter=None)
    STC.fit(data, typing['clusters'])
    
    cluster_pred = STC.transduction_
    typing['clusters'] = cluster_pred
    typing['genotype_pred'] = typing['clusters']
    del typing['clusters']
    
    return typing  # Adjust return values as needed

# Additional functions can be added here...
