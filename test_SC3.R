#### Initialisation ####
library("SingleCellExperiment")
library ("tidyverse")
library("SC3")
library("scater")
library("DescTools")
library('sva')
library("Rtsne")
library("M3C")
library("Seurat")
library("sctransform")
library("DESeq2")
library('bseqsc')

setwd("dataset/")

organise_marker_genes <- function(object, k, p_val, auroc) {
  dat <- rowData(object)[, c(paste0("sc3_", k, "_markers_clusts"), paste0("sc3_", k, 
                                                                          "_markers_auroc"), paste0("sc3_", k, "_markers_padj"), "feature_symbol")]
  dat <- dat[dat[, paste0("sc3_", k, "_markers_padj")] < p_val & !is.na(dat[, paste0("sc3_", 
                                                                                     k, "_markers_padj")]), ]
  dat <- dat[dat[, paste0("sc3_", k, "_markers_auroc")] > auroc, ]
  
  d <- NULL
  
  for (i in sort(unique(dat[, paste0("sc3_", k, "_markers_clusts")]))) {
    tmp <- dat[dat[, paste0("sc3_", k, "_markers_clusts")] == i, ]
    tmp <- tmp[order(tmp[, paste0("sc3_", k, "_markers_auroc")], decreasing = TRUE), ]
    d <- rbind(d, tmp)
  }
  
  if(nrow(dat) > 0) {
    return(d)
  } else {
    return(NULL)
  }
}

#### SC3 analysis of sce (sciatic nerves) ####

data_single <- read.table("GSE144707_countTable_aggrNerveStStD1D5.txt.gz", header = TRUE)
data_single[1:3, 1:3]
rownames(data_single) <- data_single[,1]
data_single <- subset(data_single, select = -1)

cell_class <- data.frame(cell_name = colnames(data_single))
cell_class$class <- ifelse(cell_class$cell_name %like% "%StSt", "StSt", ifelse(cell_class$cell_name %like% "%D1", "D1", "D5"))

sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_single),
    logcounts = log2(as.matrix(data_single) + 1)
  ),
  colData = cell_class
)


# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)    

plotPCA(sce, colour_by = "class")
sce <- sc3(sce, ks = 2:4, biology = TRUE)
#sce <- sc3(sce, ks = 3, biology = TRUE)

col_data <- colData(sce)
head(col_data[ , grep("sc3_", colnames(col_data))])

plotPCA(
  sce, 
  colour_by = "sc3_3_clusters", 
  size_by = "sc3_3_log2_outlier_score"
)

sc3_plot_markers(sce, k = 3)

sc3_plot_markers(
  sce, k = 3, 
  show_pdata = c(
    "class", 
    "log10_total_features",
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)



#### Lungs data analysis ####
data_lungs <- read.table("GSE113530_countsNormFinal.txt.gz", header = TRUE)

dup <- duplicated(data_lungs)
data_lungs_dedup <- data_lungs[!dup,]




#### SC3 sce2 analysis ####
sce2 <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_lungs),
    logcounts = log2(as.matrix(data_lungs) + 1)
  ))

# define feature names in feature_symbol column
rowData(sce2)$feature_symbol <- rownames(sce2)
# remove features with duplicated names
sce2 <- sce2[!duplicated(rowData(sce2)$feature_symbol), ]
# define spike-ins
isSpike(sce2, "ERCC") <- grepl("ERCC", rowData(sce2)$feature_symbol)    
plotPCA(sce2)

sce2 <- sc3(sce2, ks = NULL, biology = T, k_estimator = T)

sc3_plot_markers(
  sce2, k = 3,
  show_pdata = c(
    "sc3_3_clusters",
    "sc3_3_log2_outlier_score"
  )
)
#sc3_interactive(sce2)
marker <- organise_marker_genes(sce2, 3, 0.01, 0.85)
sc3_clust1 <- subset(marker$feature_symbol, marker$sc3_3_markers_clusts == 1)
sc3_clust3 <- subset(marker$feature_symbol,marker$sc3_3_markers_clusts == 3)
### NO MARKER FOR CLUST 2 ??????


#### sce2 unNorm ####
#test <- read.table("GSE117975_countsNormFinal.txt.gz", header = T)

data_lungs_unNorm <- read.table("GSE113530_countsFinal.txt.gz", header = TRUE)

sce2_unNorm <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_lungs_unNorm),
    logcounts = log2(as.matrix(data_lungs_unNorm) + 1)
  ))

# define feature names in feature_symbol column
rowData(sce2_unNorm)$feature_symbol <- rownames(sce2_unNorm)
# remove features with duplicated names
sce2_unNorm <- sce2_unNorm[!duplicated(rowData(sce2_unNorm)$feature_symbol), ]
# define spike-ins
isSpike(sce2_unNorm, "ERCC") <- grepl("ERCC", rowData(sce2_unNorm)$feature_symbol)    
plotPCA(sce2_unNorm)



sce2_unNorm <- sc3(sce2_unNorm, ks = NULL, biology = TRUE, k_estimator = TRUE)


sc3_interactive(sce2_unNorm)




#### SC3 analysis Tcell ####
data_Tcell <- read.table("GSE130812_FPKM_C1_table.txt.gz", header= TRUE)

rownames(data_Tcell) <- data_Tcell[,1]
gene_table <- subset(data_Tcell, select = c(2))
data_Tcell <- subset(data_Tcell, select = -c(1,2))

sce3 <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_Tcell),
    logcounts = log2(as.matrix(data_Tcell) + 1)
  ))

# define feature names in feature_symbol column
rowData(sce3)$feature_symbol <- rownames(sce3)
# remove features with duplicated names
sce3 <- sce3[!duplicated(rowData(sce3)$feature_symbol), ]
# define spike-ins
isSpike(sce3, "ERCC") <- grepl("ERCC", rowData(sce3)$feature_symbol)    
plotPCA(sce3)
sce3 <- sc3(sce3, ks = 3:6, biology = TRUE)
col_data3 <- colData(sce3)
head(col_data3[ , grep("sc3_", colnames(col_data3))])

sc3_plot_markers(sce3, k = 3)

sc3_plot_markers(
  sce2, k = 5, 
  show_pdata = c(
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)





#### SC3 analysis broad digestive organs ####
data_broad <- read.table("GSE92332_AtlasFullLength_TPM.txt.gz",header = TRUE)

sce4 <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data_broad),
    logcounts = log2(as.matrix(data_broad) + 1)
  ))

rowData(sce4)$feature_symbol <- rownames(sce4)
# remove features with duplicated names
sce4 <- sce4[!duplicated(rowData(sce4)$feature_symbol), ]
# define spike-ins
isSpike(sce4, "ERCC") <- grepl("ERCC", rowData(sce4)$feature_symbol)    
plotPCA(sce4)

sce4 <- sc3(sce4, ks = 3:6, biology = TRUE)

col_data4 <- colData(sce4)
head(col_data4[ , grep("sc3_", colnames(col_data4))])
plotPCA(sce4, colour_by = "sc3_6_clusters")

sc3_plot_markers(sce4, k = 3)

sc3_plot_markers(
  sce4, k = 3, 
  show_pdata = c(
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)


marker <- organise_marker_genes(sce4, 4, 0.01, 0.80)

marker_clust1 <- get_marker_genes(sce4, )


#### Seurat lungs ####

seurat_lungs_data <- CreateSeuratObject(counts = data_lungs_unNorm)
seurat_lungs_data <- SCTransform(seurat_lungs_data)

seurat_lungs_data <- RunPCA(seurat_lungs_data)
seurat_lungs_data <- RunUMAP(seurat_lungs_data, dims = 1:30)
seurat_lungs_data <- FindNeighbors(seurat_lungs_data, dims = 1:30)
seurat_lungs_data <- FindClusters(seurat_lungs_data)
DimPlot(seurat_lungs_data, label = T) + NoLegend()
seurat_markers <- FindAllMarkers(seurat_lungs_data, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)
seurat_final_markers <- seurat_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
seurat_clust1 <- subset(seurat_final_markers$gene, seurat_final_markers$cluster == 1)
seurat_clust0 <- subset(seurat_final_markers$gene, seurat_final_markers$cluster == 0)


seurat_lungs_Norm_data <- CreateSeuratObject(counts = data_lungs)
seurat_lungs_Norm_data <- SCTransform(seurat_lungs_Norm_data)

seurat_lungs_Norm_data <- RunPCA(seurat_lungs_Norm_data)
seurat_lungs_Norm_data <- RunUMAP(seurat_lungs_Norm_data, dims = 1:30)
seurat_lungs_Norm_data <- FindNeighbors(seurat_lungs_Norm_data, dims = 1:30)
seurat_lungs_Norm_data <- FindClusters(seurat_lungs_Norm_data)
DimPlot(seurat_lungs_Norm_data, label = T) + NoLegend()
seurat_markers_Norm <- FindAllMarkers(seurat_lungs_Norm_data, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)
seurat_final_markers_Norm <- seurat_markers_Norm %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
seurat_clust1_Norm <- subset(seurat_final_markers_Norm$gene, seurat_final_markers_Norm$cluster == 1)
seurat_clust0_Norm <- subset(seurat_final_markers_Norm$gene, seurat_final_markers_Norm$cluster == 0)

#### Seurat Sciatic Nerves ####

siatic <- CreateSeuratObject(counts = data_single)
siatic <- SCTransform(siatic)
siatic <- RunPCA(siatic)
siatic <- RunUMAP(siatic, dims = 1:30)
siatic <- FindNeighbors(siatic, dims = 1:30)
siatic <- FindClusters(siatic)
DimPlot(siatic, label = T) + NoLegend()
siatic_markers <- FindAllMarkers(siatic, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)
siatic_markers_final <- siatic_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

write_delim(siatic_markers_final, "seurat_markers_GSE144707", "\t")



#### Seurat broad####
seurat_broad_data <- CreateSeuratObject(counts = data_broad)
seurat_broad_data <- SCTransform(seurat_broad_data)

seurat_broad_data <- RunPCA(seurat_broad_data)
seurat_broad_data <- RunUMAP(seurat_broad_data, dims = 1:30)
seurat_broad_data <- FindNeighbors(seurat_broad_data, dims = 1:30)
seurat_broad_data <- FindClusters(seurat_broad_data)
DimPlot(seurat_broad_data, label = T) + NoLegend()
seurat_markers_broad <- FindAllMarkers(seurat_broad_data, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)
seurat_final_markers_broad <- seurat_markers_broad %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

write_delim(seurat_final_markers_broad, "seurat_markers_GSE92332", "\t")
#### comparaison marker ####

seurat_clust0 <- gsub("-","_", seurat_clust0)
seurat_clust0_Norm <- gsub("-","_", seurat_clust0_Norm)
seurat_clust1 <- gsub("-","_", seurat_clust1)
seurat_clust1_Norm <- gsub("-","_", seurat_clust1_Norm)
intersect(seurat_clust0,sc3_clust1)
intersect(seurat_clust0,sc3_clust3)
intersect(seurat_clust0, seurat_clust0_Norm)
intersect(seurat_clust0, seurat_clust1_Norm)
intersect(seurat_clust1,sc3_clust1)
intersect(seurat_clust1,sc3_clust3)
intersect(seurat_clust1, seurat_clust0_Norm)
intersect(seurat_clust1, seurat_clust1_Norm)

intersect(seurat_clust0_Norm, sc3_clust1)
intersect(seurat_clust0_Norm, sc3_clust3)
intersect(seurat_clust1_Norm, sc3_clust1)
intersect(seurat_clust1_Norm, sc3_clust3)


write_delim(as.data.frame(seurat_clust0), 'seurat_clust0_markers_GSE113530', '\t')
write_delim(as.data.frame(seurat_clust1), 'seurat_clust1_markers_GSE113530', '\t')


#### Bulk sciatic ####

data_bulk <- read.table("GSE144705_processedData_bulkRNAseq_YdensEtAl.txt.gz", header = TRUE)
data_bulk

data_bulk_marker <- subset(data_bulk, data_bulk$Gene %like% siatic_markers$gene)
rownames(data_bulk_marker) <- data_bulk_marker$Gene
data_bulk_marker <- subset(data_bulk_marker, select= -1)
bulk_samples <- data.frame(matrix(nrow = 13 ,ncol =0 ))
rownames(bulk_samples) <- colnames(data_bulk_marker)[2:14]
bulk_samples$condition <- rep(c("SN","ON","SPF"),c(4,4,5))

dds <- DESeqDataSetFromMatrix(countData = data_bulk_marker,
                              colData = bulk_samples,
                              design = ~ condition)


dds <- DESeq(dds)
res_SN_ON <- results(dds, contrast = c("condition","SN","ON"))
res_SN_SPF <- results(dds, contrast = c("condition","SN","SPF"))
res_ON_SPF <- results(dds, contrast = c("condition","SPF","ON"))

sum(res$padj < 0.1, na.rm=TRUE)
sum(res_ON_SPF$padj < 0.1, na.rm=TRUE)
sum(res_SN_ON$padj < 0.1, na.rm=TRUE)
sum(res_SN_SPF$padj < 0.1, na.rm=TRUE)



#### Pancreas dataset ####
pancreas_h1 <- read.csv("GSE84133_Pancreas/GSM2230757_human1_umifm_counts.csv.gz", header = TRUE)
rownames(pancreas_h1) <- pancreas_h1$X
pancreas_h1 <- subset(pancreas_h1, select = -c(1))
pancreas_h1_sup <- subset(pancreas_h1, select = c(1,2))
pancreas_h1 <- subset(pancreas_h1, select = -c(1,2))
pancreas_h1_long <- pancreas_h1
pancreas_h1 <- t(pancreas_h1)

pancreas_h1_seurat <- CreateSeuratObject(counts = pancreas_h1)
pancreas_h1_seurat <- SCTransform(pancreas_h1_seurat)

pancreas_h1_seurat <- RunPCA(pancreas_h1_seurat)
pancreas_h1_seurat <- RunUMAP(pancreas_h1_seurat, dims = 1:30)
pancreas_h1_seurat <- FindNeighbors(pancreas_h1_seurat, dims = 1:30)
pancreas_h1_seurat <- FindClusters(pancreas_h1_seurat)
DimPlot(pancreas_h1_seurat, label = T) + NoLegend()
pancreas_h1_seurat_markers <- FindAllMarkers(pancreas_h1_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)

pancreas_h1_seurat <- RunTSNE(pancreas_h1_seurat)
TSNEPlot(pancreas_h1_seurat)


pancreas_h2 <- read.csv("GSE84133_Pancreas/GSM2230758_human2_umifm_counts.csv.gz", header = TRUE)
rownames(pancreas_h2) <- pancreas_h2$X
pancreas_h2 <- subset(pancreas_h2, select = -c(1))
pancreas_h2_sup <- subset(pancreas_h2, select = c(1,2))
pancreas_h2 <- subset(pancreas_h2, select = -c(1,2))
pancreas_h2_long <- pancreas_h2
pancreas_h2 <- t(pancreas_h2)

pancreas_h2_seurat <- CreateSeuratObject(counts = pancreas_h2)
pancreas_h2_seurat <- SCTransform(pancreas_h2_seurat)

pancreas_h2_seurat <- RunPCA(pancreas_h2_seurat)
pancreas_h2_seurat <- RunUMAP(pancreas_h2_seurat, dims = 1:30)
pancreas_h2_seurat <- FindNeighbors(pancreas_h2_seurat, dims = 1:30)
pancreas_h2_seurat <- FindClusters(pancreas_h2_seurat)
DimPlot(pancreas_h2_seurat, label = T) + NoLegend()
pancreas_h2_seurat_markers <- FindAllMarkers(pancreas_h2_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)





pancreas_h3 <- read.csv("GSE84133_Pancreas/GSM2230759_human3_umifm_counts.csv.gz", header = TRUE)
rownames(pancreas_h3) <- pancreas_h3$X
pancreas_h3 <- subset(pancreas_h3, select = -c(1))
pancreas_h3_sup <- subset(pancreas_h3, select = c(1,2))
pancreas_h3 <- subset(pancreas_h3, select = -c(1,2))
pancreas_h3_long <- pancreas_h3
pancreas_h3 <- t(pancreas_h3)

pancreas_h3_seurat <- CreateSeuratObject(counts = pancreas_h3)
pancreas_h3_seurat <- SCTransform(pancreas_h3_seurat)

pancreas_h3_seurat <- RunPCA(pancreas_h3_seurat)
pancreas_h3_seurat <- RunUMAP(pancreas_h3_seurat, dims = 1:30)
pancreas_h3_seurat <- FindNeighbors(pancreas_h3_seurat, dims = 1:30)
pancreas_h3_seurat <- FindClusters(pancreas_h3_seurat)
DimPlot(pancreas_h3_seurat, label = T) + NoLegend()
pancreas_h3_seurat_markers <- FindAllMarkers(pancreas_h3_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)




pancreas_h4 <- read.csv("GSE84133_Pancreas/GSM2230760_human4_umifm_counts.csv.gz", header = TRUE)
rownames(pancreas_h4) <- pancreas_h4$X
pancreas_h4 <- subset(pancreas_h4, select = -c(1))
pancreas_h4_sup <- subset(pancreas_h4, select = c(1,2))
pancreas_h4 <- subset(pancreas_h4, select = -c(1,2))
pancreas_h4_long <- pancreas_h4
pancreas_h4 <- t(pancreas_h4)

pancreas_h4_seurat <- CreateSeuratObject(counts = pancreas_h4)
pancreas_h4_seurat <- SCTransform(pancreas_h4_seurat)

pancreas_h4_seurat <- RunPCA(pancreas_h4_seurat)
pancreas_h4_seurat <- RunUMAP(pancreas_h4_seurat, dims = 1:30)
pancreas_h4_seurat <- FindNeighbors(pancreas_h4_seurat, dims = 1:30)
pancreas_h4_seurat <- FindClusters(pancreas_h4_seurat)
DimPlot(pancreas_h4_seurat, label = T) + NoLegend()
pancreas_h4_seurat_markers <- FindAllMarkers(pancreas_h4_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)





pancreas_all <-rbind(pancreas_h1_long,pancreas_h2_long,pancreas_h3_long,pancreas_h4_long)
pancreas_all <- t(pancreas_all)

pancreas_all_seurat <- CreateSeuratObject(counts = pancreas_all)
pancreas_all_seurat <- SCTransform(pancreas_all_seurat)

pancreas_all_seurat <- RunPCA(pancreas_all_seurat)
pancreas_all_seurat <- RunTSNE(pancreas_all_seurat)
TSNEPlot(pancreas_all_seurat)
pancreas_all_seurat <- RunUMAP(pancreas_all_seurat, dims = 1:30)
pancreas_all_seurat <- FindNeighbors(pancreas_all_seurat, dims = 1:30)
pancreas_all_seurat <- FindClusters(pancreas_all_seurat)
DimPlot(pancreas_all_seurat, label = T) + NoLegend()
pancreas_all_seurat_markers <- FindAllMarkers(pancreas_all_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)

pancreas_all_seurat <- JackStraw(pancreas_all_seurat, num.replicate = 100)
pancreas_all_seurat <-  ScoreJackStraw(pancreas_all_seurat, dims = 1:20)

JackStrawPlot(pancreas_all_seurat, dims = 1:20)

ElbowPlot(pancreas_all_seurat,)



length(intersect(pancreas_h2_seurat_markers$gene, pancreas_h1_seurat_markers$gene))
length(pancreas_h1_seurat_markers$gene)
length(pancreas_h2_seurat_markers$gene)

known_marker <- c("GCG",'INS','PPY','SST','GHRL','PRSS1','KRT19','SPARC','VWF','RGS5','PDGFRA','SOX10','SDS','TPSAB1','TRAC')

topn_all <- pancreas_all_seurat_markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
DoHeatmap(pancreas_all_seurat, features = known_marker) + NoLegend()

topn_h1 <- pancreas_h1_seurat_markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
DoHeatmap(pancreas_h1_seurat, features = known_marker) + NoLegend()
DoHeatmap(pancreas_h2_seurat, features = known_marker) + NoLegend()
DoHeatmap(pancreas_h3_seurat, features = known_marker) + NoLegend()
DoHeatmap(pancreas_h4_seurat, features = known_marker) + NoLegend()





