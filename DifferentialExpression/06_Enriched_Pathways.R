#################### ROntoTools ############################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ROntoTools")
# 
# 
# require(graph)
# require(ROntoTools)
# # Get pathway graphs and names
# kpg <- keggPathwayGraphs("cge", updateCache = TRUE)
# kpn <- keggPathwayNames("cge")
# 
# # Analyze the first sample
# wtVS6xDay4 <- read.csv("Line1_SignificantGenes.csv")
# View(head(wtVS6xDay4))
# 
# # Prepare vectors for analysis
# fc <- wtVS6xDay4$log2FoldChange
# names(fc) <- wtVS6xDay4$keggGeneName
# pv <- wtVS6xDay4$padj
# names(pv) <- wtVS6xDay4$keggGeneName
# ref <- wtVS6xDay4$keggGeneName
# 
# # Set node weights in KEGG graph
# kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
# 
# # Run the Pathway Enrichment analysis
# peRes <- pe(x = fc, graphs = kpg, ref = ref, nboot = 500, verbose = TRUE)
# 
# # Summarize results
# Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE, pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert")
# 
# # Plot results
# plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = .01)
# plot(peRes@pathways[["path:cge04072"]], type = "two.way")
# plot(peRes@pathways[["path:cge04072"]], type = "boot")
# 
# # Render graphs
# p <- peRes@pathways[["path:cge04072"]]
# g <- layoutGraph(p@map, layoutType = "dot")
# graphRenderInfo(g) <- list(fixedsize = FALSE)
# edgeRenderInfo(g) <- peEdgeRenderInfo(p)
# nodeRenderInfo(g) <- peNodeRenderInfo(p)
# renderGraph(g)

################## Fast GSEA #######################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("fgsea")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("reactome.db")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")
setwd("C:\\Users\\Jahir\ Gutierrez\\Documents\\R\\DiffExp_StefPaper")
library(fgsea)
library(graph)
library(ROntoTools)
library(reactome.db)
library(biomaRt)
kpg <- keggPathwayGraphs("cge", updateCache = TRUE)
kpn <- keggPathwayNames("cge")
data("examplePathways")
data("exampleRanks")

