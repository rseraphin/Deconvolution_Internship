#### Initialisation ####
if (!require('tidyverse')){
  install.packages("tidyverse")
  library ("tidyverse")
}
if (!require('Seurat')){
  BiocManager::install("multtest")
  install.packages('Seurat')
  library("Seurat")
  install.packages('sctransform')
  library("sctransform")
}
if (!require('nnls')){
  install.packages('nnls')
  library('nnls')
}

setwd("Master_2/")
#load('dataset/.RData')

## Reading and cleaning the data
pancreas_h1 <- read.csv("dataset/GSE84133_Pancreas/GSM2230757_human1_umifm_counts.csv.gz", header = TRUE)
rownames(pancreas_h1) <- pancreas_h1$X
pancreas_h1 <- subset(pancreas_h1, select = -c(1))
pancreas_h1_sup <- subset(pancreas_h1, select = c(1,2))
pancreas_h1 <- subset(pancreas_h1, select = -c(1,2))
pancreas_h1_long <- pancreas_h1
pancreas_h1 <- t(pancreas_h1)

## Creating Seurat Object
pancreas_h1_seurat <- CreateSeuratObject(counts = pancreas_h1)
## Performing Normalization, Imputation, Scaling
pancreas_h1_seurat <- SCTransform(pancreas_h1_seurat)

## Studiying dimensionality
pancreas_h1_seurat <- RunPCA(pancreas_h1_seurat)
pancreas_h1_seurat <- JackStraw(pancreas_h1_seurat, num.replicate = 100, dims = 30)
pancreas_h1_seurat <-  ScoreJackStraw(pancreas_h1_seurat, dims = 1:30)

JackStrawPlot(pancreas_h1_seurat, dims = 1:30)
ElbowPlot(pancreas_h1_seurat,ndims = 30)

## Runnig dimension reduction techniques
pancreas_h1_seurat <- RunUMAP(pancreas_h1_seurat, dims = 1:30)

## Clustering
pancreas_h1_seurat <- FindNeighbors(pancreas_h1_seurat, dims = 1:30)
#pancreas_h1_seurat <- FindClusters(pancreas_h1_seurat, resolution = 0.6)
pancreas_h1_seurat <- FindClusters(pancreas_h1_seurat)
DimPlot(pancreas_h1_seurat, label = T) + NoLegend()

## Finding markers
pancreas_h1_seurat_markers <- FindAllMarkers(pancreas_h1_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)


## Exporting data files
# Meta file
meta_data_h1 <- pancreas_h1_seurat@meta.data

# Raw sc data
treated_h1_data <- as.matrix(GetAssayData(pancreas_h1_seurat,slot='data'))

# Avg sc data
pancreas_h1_avg <- AverageExpression(pancreas_h1_seurat)

write.table(meta_data_h1, 'dataset/meta_data_h1', quote = FALSE,sep = '\t', col.names = NA)
write.table(as.data.frame(treated_h1_data), 'dataset/treated_h1_data', quote = FALSE, sep = '\t', col.names = NA)

## Data for deconvolution
# bulk data
pancreas_bulk <- read.table('dataset/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz', header = T)

selected_genes_bulk_h1 <- subset(pancreas_bulk, pancreas_bulk$id %in% pancreas_h1_seurat_markers$gene)
selected_genes_sc_h1 <- subset(pancreas_h1_avg$SCT, rownames(pancreas_h1_avg$SCT) %in% selected_genes_bulk_h1$id)

write_tsv(selected_genes_bulk_h1, 'dataset/selected_genes_bulk_h1')
write.table(selected_genes_sc_h1, 'dataset/selected_genes_sc_h1', sep='\t' , quote=FALSE, col.names = NA)
