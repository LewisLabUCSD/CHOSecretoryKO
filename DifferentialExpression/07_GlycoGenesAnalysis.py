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

f = open('HyperGeometric_Test_Results_GlycoGenes.txt', 'w')
# Run HyperGeometric tests
for i in range(len(DE_files)):
    description = DE_desc[DE_files[i]]
    print("Processing file {} - {}".format(DE_files[i],description))
    glyco_genes = pd.read_csv("GlycoGenes_CHO.csv", header = 0, dtype={'geneid':str})
    entrezToName = pd.read_csv("EntrezToNameMap.csv", header=0)
    diff_expression = pd.read_csv(DE_files[i], header = 0)
    diff_expression = diff_expression.rename(columns={"Unnamed: 0":"gename"})
    diff_expression = pd.merge(entrezToName, diff_expression, on="gename")
    common_genes = list(set(glyco_genes['geneid']) & set(diff_expression['geneid']))
    idx = [i for i in range(len(glyco_genes)) if glyco_genes['geneid'][i] in common_genes]
    glyco_genes = glyco_genes.iloc[idx,:].reset_index(drop=True)
    idx = [i for i in range(len(diff_expression)) if isSignificant(diff_expression['log2FoldChange'][i], diff_expression['padj'][i])]
    significant_genes = diff_expression.iloc[idx,:].reset_index(drop=True)
    P = len(set(diff_expression['geneid'])) # Population
    S = len(set(significant_genes['geneid'])) # Successes in population
    p = len(set(glyco_genes['geneid'])) # Sample
    s = len(set(glyco_genes['geneid']) & set(significant_genes['geneid'])) # Successess in sample
    f.write('SAMPLE COMPARISON (All glycosyltransferase genes): {}\n'.format(description))
    f.write('Total genes: {}\n'.format(P))
    f.write('Differentially expressed genes: {}\n'.format(S))
    f.write('Number of glyco genes: {}\n'.format(p))
    f.write('Differentially expressed glyco genes (s): {}\n'.format(s))
    f.write('p-value s <= ' + str(s) + ': ' + str(hypergeom.cdf(s ,P,S,p)) + '\n')
    f.write('p-value s >= ' + str(s) + ': ' + str(hypergeom.sf(s - 1,P,S,p)) + '\n\n')
    for group in list(set(glyco_genes['class'])):
        P = len(set(diff_expression['geneid'])) # Population
        S = len(set(significant_genes['geneid'])) # Successes in population
        p = len(set(glyco_genes[glyco_genes['class']==group]['geneid'])) # Sample
        s = len(set(glyco_genes[glyco_genes['class']==group]['geneid']) & set(significant_genes['geneid'])) # Successess in sample
        f.write('SAMPLE COMPARISON ({}): {}\n'.format(group,description))
        f.write('Total genes: {}\n'.format(P))
        f.write('Differentially expressed genes: {}\n'.format(S))
        f.write('Number of glyco genes: {}\n'.format(p))
        f.write('Differentially expressed glyco genes (s): {}\n'.format(s))
        f.write('p-value s <= ' + str(s) + ': ' + str(hypergeom.cdf(s ,P,S,p)) + '\n')
        f.write('p-value s >= ' + str(s) + ': ' + str(hypergeom.sf(s - 1,P,S,p)) + '\n\n')
f.close()

print("Finished with success!")