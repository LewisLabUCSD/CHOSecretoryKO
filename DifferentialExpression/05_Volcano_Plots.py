import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

csv_files = os.listdir(os.getcwd())
csv_files = [f for f in csv_files if "Line" in f and ".csv" in f]

# Function to determine significance
def isSignificant(xval,yval, xthr = 1, ythr = 2):
    if abs(xval) >= xthr and abs(yval) >= ythr:
        return True
    else:
        return False

# Read Entrez -> Name map
entrezToName = pd.read_csv("EntrezToNameMap.csv", header=0)

for csv_file in csv_files:
    print("Processing file {}".format(csv_file))
    df = pd.read_csv(csv_file, header=0)
    df = df.rename(columns={"Unnamed: 0":"gename"})
    x = df['log2FoldChange'].values
    y = df['padj'].values + 1e-5
    y = -np.log10(y)
    significant_idx = [i for i in range(len(x)) if isSignificant(x[i],y[i])]
    nonsignificant_idx = [i for i in range(len(x)) if not isSignificant(x[i],y[i])]

    # Plot Volcano Plot
    plt.figure(figsize=(8,8))
    plt.scatter(x[significant_idx], y[significant_idx], c='red', alpha=0.35, label='Significant')
    plt.scatter(x[nonsignificant_idx], y[nonsignificant_idx], c='blue', alpha=0.35, label='Nonsignificant')
    plt.vlines(-1, 0, 5, linestyles='dashed')
    plt.vlines(1, 0, 5, linestyles='dashed')
    plt.hlines(2, min(x), max(x), linestyles='dashed')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10 (adjusted p-value)')
    plt.legend()
    plt.savefig(csv_file.replace(".csv","_volcanoPlot.pdf"))

    # Save names of significant differentially expressed genes
    tmp_df = df.iloc[significant_idx,:].reset_index(drop=True)
    final_df = pd.merge(entrezToName, tmp_df, on="gename")
    final_df['keggGeneName'] = ["cge:" + str(id) for id in list(final_df['geneid'])] # Required for pathway analysis with ROntoTools
    final_df.to_csv(csv_file.replace(".csv","_SignificantGenes.csv"), index=False)