import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom

# Function to assess significant differential expression
def isSignificant(log2FC, padj):
    if abs(log2FC) >= 1 and -np.log10(padj + 1e-5) >= 2:
        return True
    else:
        return False

# Prepare list of files
DE_files = ['Line{}.csv'.format(i) for i in range(1,7)]
DE_desc = {'Line1.csv':'WT vs 6X day 4', 'Line2.csv':'WT vs 6X day 6', 'Line3.csv':'WT vs 6X day 8', 'Line4.csv':'WT vs 11X day 4', 'Line5.csv':'WT vs 11X day 6', 'Line6.csv':'WT vs 11X day 8'}

f = open('HyperGeometric_Test_Results_SecretoryGenes.txt', 'w')
# Run HyperGeometric tests
for i in range(len(DE_files)):
    description = DE_desc[DE_files[i]]
    print("Processing file {} - {}".format(DE_files[i],description))
    secretory_genes = pd.read_csv("SecretoryGenes_CHO.csv", header = 0, dtype={'geneid':str})
    entrezToName = pd.read_csv("EntrezToNameMap.csv", header=0)
    diff_expression = pd.read_csv(DE_files[i], header = 0)
    diff_expression = diff_expression.rename(columns={"Unnamed: 0":"gename"})
    diff_expression = pd.merge(entrezToName, diff_expression, on="gename")
    common_genes = list(set(secretory_genes['geneid']) & set(diff_expression['geneid']))
    idx = [i for i in range(len(secretory_genes)) if secretory_genes['geneid'][i] in common_genes]
    secretory_genes = secretory_genes.iloc[idx,:].reset_index(drop=True)
    idx = [i for i in range(len(diff_expression)) if isSignificant(diff_expression['log2FoldChange'][i], diff_expression['padj'][i])]
    significant_genes = diff_expression.iloc[idx,:].reset_index(drop=True)
    P = len(set(diff_expression['geneid'])) # Population
    S = len(set(significant_genes['geneid'])) # Successes in population
    p = len(set(secretory_genes['geneid'])) # Sample
    s = len(set(secretory_genes['geneid']) & set(significant_genes['geneid'])) # Successess in sample
    f.write('SAMPLE COMPARISON (All secretory genes): {}\n'.format(description))
    f.write('Total genes: {}\n'.format(P))
    f.write('Differentially expressed genes: {}\n'.format(S))
    f.write('Number of secretory genes: {}\n'.format(p))
    f.write('Differentially expressed secretory genes (s): {}\n'.format(s))
    f.write('p-value s <= ' + str(s) + ': ' + str(hypergeom.cdf(s ,P,S,p)) + '\n')
    f.write('p-value s >= ' + str(s) + ': ' + str(hypergeom.sf(s - 1,P,S,p)) + '\n\n')
    for group in list(set(secretory_genes['module'])):
        P = len(set(diff_expression['geneid'])) # Population
        S = len(set(significant_genes['geneid'])) # Successes in population
        p = len(set(secretory_genes[secretory_genes['module']==group]['geneid'])) # Sample
        s = len(set(secretory_genes[secretory_genes['module']==group]['geneid']) & set(significant_genes['geneid'])) # Successess in sample
        pval_left = hypergeom.cdf(s ,P,S,p)
        pval_right = hypergeom.sf(s - 1,P,S,p)
        if pval_left <= 0.05 or pval_right <= 0.05:
            f.write('SAMPLE COMPARISON ({}): {}\n'.format(group,description))
            f.write('Total genes: {}\n'.format(P))
            f.write('Differentially expressed genes: {}\n'.format(S))
            f.write('Number of secretory genes: {}\n'.format(p))
            f.write('Differentially expressed secretory genes (s): {}\n'.format(s))
            f.write('p-value s <= ' + str(s) + ': ' + str(hypergeom.cdf(s ,P,S,p)) + '\n')
            f.write('p-value s > ' + str(s) + ': ' + str(hypergeom.sf(s - 1,P,S,p)) + '\n\n')
        else:
            continue
f.close()

print("Finished with success!")