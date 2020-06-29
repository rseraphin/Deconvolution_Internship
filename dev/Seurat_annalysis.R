
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

#### Seurat lungs ####
data_lungs_unNorm <- read.table("dataset/GSE113530_countsFinal.txt", header = TRUE)
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
data_single <- read.table("dataset/GSE144707_countTable_aggrNerveStStD1D5.txt", header = TRUE)
data_single[1:3, 1:3]
rownames(data_single) <- data_single[,1]
data_single <- subset(data_single, select = -1)


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

data_bulk <- read.table("dataset/GSE144705_processedData_bulkRNAseq_YdensEtAl.txt.gz", header = TRUE)
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
pancreas_h1 <- read.csv("dataset/GSE84133_Pancreas/GSM2230757_human1_umifm_counts.csv.gz", header = TRUE)
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
#pancreas_h1_seurat <- FindClusters(pancreas_h1_seurat, resolution = 0.6)
pancreas_h1_seurat <- FindClusters(pancreas_h1_seurat)
DimPlot(pancreas_h1_seurat, label = T) + NoLegend()
pancreas_h1_seurat_markers <- FindAllMarkers(pancreas_h1_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)
pancreas_h1_seurat_markers_c <- FindAllMarkers(pancreas_h1_seurat, only.pos = T)

pancreas_h1_seurat <- RunTSNE(pancreas_h1_seurat, )
TSNEPlot(pancreas_h1_seurat)

pancreas_h1_seurat <- JackStraw(pancreas_h1_seurat, num.replicate = 100, dims = 20)
pancreas_h1_seurat <-  ScoreJackStraw(pancreas_h1_seurat, dims = 1:20)

JackStrawPlot(pancreas_h1_seurat, dims = 1:20)

ElbowPlot(pancreas_h1_seurat,ndims = 20)


meta_data_h1 <- pancreas_h1_seurat@meta.data

treated_h1_data <- as.matrix(GetAssayData(pancreas_h1_seurat,slot='data'))

pancreas_h1_avg <- AverageExpression(pancreas_h1_seurat)

write.table(meta_data_h1, 'dataset/meta_data_h1', quote = FALSE,sep = '\t', col.names = NA)
write.table(as.data.frame(treated_h1_data), 'dataset/treated_h1_data', quote = FALSE, sep = '\t', col.names = NA)


selected_genes_bulk_h1 <- subset(pancreas_bulk, pancreas_bulk$id %in% pancreas_h1_seurat_markers$gene)
selected_genes_sc_h1 <- subset(pancreas_h1_avg$SCT, rownames(pancreas_h1_avg$SCT) %in% selected_genes_bulk_h1$id)

write_tsv(selected_genes_bulk_h1, 'dataset/selected_genes_bulk_h1')
write.table(selected_genes_sc_h1, 'dataset/selected_genes_sc_h1', sep='\t' , quote=FALSE, col.names = NA)

pancreas_h2 <- read.csv("dataset/GSE84133_Pancreas/GSM2230758_human2_umifm_counts.csv.gz", header = TRUE)
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





meta_data_h2 <- pancreas_h2_seurat@meta.data

treated_h2_data <- as.matrix(GetAssayData(pancreas_h2_seurat,slot='data'))

pancreas_h2_avg <- AverageExpression(pancreas_h2_seurat)

write.table(meta_data_h2, 'dataset/meta_data_h2', quote = FALSE,sep = '\t', col.names = NA)
write.table(as.data.frame(treated_h2_data), 'dataset/treated_h2_data', quote = FALSE, sep = '\t', col.names = NA)


selected_genes_bulk_h2 <- subset(pancreas_bulk, pancreas_bulk$id %in% pancreas_h2_seurat_markers$gene)
selected_genes_sc_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% selected_genes_bulk_h2$id)

write_tsv(selected_genes_bulk_h2, 'dataset/selected_genes_bulk_h2')
write.table(selected_genes_sc_h2, 'dataset/selected_genes_sc_h2', sep='\t' , quote=FALSE, col.names = NA)







pancreas_h3 <- read.csv("dataset/GSE84133_Pancreas/GSM2230759_human3_umifm_counts.csv.gz", header = TRUE)
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




pancreas_h4 <- read.csv("dataset/GSE84133_Pancreas/GSM2230760_human4_umifm_counts.csv.gz", header = TRUE)
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
pancreas_all_seurat <- FindClusters(pancreas_all_seurat, resolution = 1.5)
DimPlot(pancreas_all_seurat, label = T) + NoLegend()
pancreas_all_seurat_markers <- FindAllMarkers(pancreas_all_seurat, only.pos = T, min.pct =  0.25, logfc.threshold = 0.25)


DoHeatmap(pancreas_all_seurat, features = known_marker) + NoLegend()

# Tested resolution 0.5 ~ 1.3

# Testing merging clusters


pancreas_all_seurat <- BuildClusterTree(pancreas_all_seurat, do.reorder = T, reorder.numeric = T)
node.scores <- AssessNodes(pancreas_all_seurat)
node.scores[order(node.scores$oobe,decreasing = T),] -> node.scores



# Components analysis
pancreas_all_seurat <- JackStraw(pancreas_all_seurat, num.replicate = 100, dims = 30)
pancreas_all_seurat <-  ScoreJackStraw(pancreas_all_seurat, dims = 1:30)

JackStrawPlot(pancreas_all_seurat, dims = 10:30)

ElbowPlot(pancreas_all_seurat,ndims = 30)


# Markers comparison and visualization
length(intersect(pancreas_h2_seurat_markers$gene, pancreas_h1_seurat_markers$gene))
length(pancreas_h1_seurat_markers$gene)
length(pancreas_h2_seurat_markers$gene)

known_marker <- c("GCG",'INS','PPY','SST','GHRL','PRSS1',"CPA1",'KRT19','SPARC','VWF','RGS5','PDGFRA','SOX10','SDS','TPSAB1','TRAC')
cell_types <- c('Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Acinar', 'Acinar', 'Ductal', '', 'Endothelial', 'quiescent stellate', 'activate stellate', 'schwann', 'Macrophage', 'Mast', 'Cytotoxyc T')
topn_all <- pancreas_all_seurat_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
DoHeatmap(pancreas_all_seurat, features = known_marker) + NoLegend()

known_marker %in% pancreas_all_seurat_markers$gene
topn_all$gene %in% pancreas_all_seurat

subset(pancreas_h1_seurat_markers, pancreas_h1_seurat_markers$gene %in% known_marker) 
subset(pancreas_all_seurat_markers, pancreas_all_seurat_markers$gene %in% known_marker) 



topn_h1 <- pancreas_h1_seurat_markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
top_1_h2 <- pancreas_h2_seurat_markers %>% group_by(cluster) %>% top_n(n=1, wt = avg_logFC)
top_5_h2 <- pancreas_h2_seurat_markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)
top_10_h2 <- pancreas_h2_seurat_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_logFC)
top_20_h2 <- pancreas_h2_seurat_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_logFC)
top_50_h2 <- pancreas_h2_seurat_markers %>% group_by(cluster) %>% top_n(n=50, wt = avg_logFC)
DefaultAssay(object = pancreas_h1_seurat) <- "SCT"
t1 <- DoHeatmap(pancreas_h1_seurat, features = known_marker) + NoLegend()
t2 <- DoHeatmap(pancreas_h2_seurat, features = known_marker) + NoLegend()
t3 <- DoHeatmap(pancreas_h3_seurat, features = known_marker) + NoLegend()
t4 <- DoHeatmap(pancreas_h4_seurat, features = known_marker) + NoLegend()

par(mfrow=c(2,2))
plot(t2)
DoHeatmap(pancreas_all_seurat, features = topn_all) + NoLegend()

top1_h2_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% top_1_h2$gene)
top5_h2_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% top_5_h2$gene)
top10_h2_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% top_10_h2$gene)
top20_h2_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% top_20_h2$gene)
top50_h2_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% top_50_h2$gene)

write_tsv(top1_h2_gene_bulk, "dataset/top1_h2_genes_bulk.tsv")
write_tsv(top5_h2_gene_bulk, "dataset/top5_h2_genes_bulk.tsv")
write_tsv(top10_h2_gene_bulk, "dataset/top10_h2_genes_bulk.tsv")
write_tsv(top20_h2_gene_bulk, "dataset/top20_h2_genes_bulk.tsv")
write_tsv(top50_h2_gene_bulk, "dataset/top50_h2_genes_bulk.tsv")

top1_genes_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% top1_h2_gene_bulk$id)
top5_genes_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% top5_h2_gene_bulk$id)
top10_genes_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% top10_h2_gene_bulk$id)
top20_genes_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% top20_h2_gene_bulk$id)
top50_genes_h2 <- subset(pancreas_h2_avg$SCT, rownames(pancreas_h2_avg$SCT) %in% top50_h2_gene_bulk$id)

write_tsv(top1_genes_h2, "dataset/top1_genes_h2.tsv")
write_tsv(top10_genes_h2, "dataset/top10_genes_h2.tsv")
write_tsv(top5_genes_h2, "dataset/top5_genes_h2.tsv")
write_tsv(top20_genes_h2, "dataset/top20_genes_h2.tsv")
write_tsv(top50_genes_h2, "dataset/top50_genes_h2.tsv")


pancreas_avg <- AverageExpression(pancreas_all_seurat)


#selected_genes_sc <- subset(pancreas_avg$SCT, rownames(pancreas_avg$SCT) %in% pancreas_all_seurat_markers$gene)

GetAssayData(object = pancreas_all_seurat, slot = 'data')[1:3,1:3]
GetAssayData(object = pancreas_h1_seurat, slot = 'RNA')
#### Bulk Pancreas GSE50244 ####

pancreas_bulk <- read.table('dataset/GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz', header = T)

selected_genes_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% pancreas_all_seurat_markers$gene)
selected_genes_sc <- subset(pancreas_avg$SCT, rownames(pancreas_avg$SCT) %in% selected_genes_bulk$id)

top20_gene_bulk <- subset(pancreas_bulk, pancreas_bulk$id %in% topn_all$gene)
top20_genes_sc <- subset(pancreas_avg$SCT, rownames(pancreas_avg$SCT) %in% top20_gene_bulk$id)

write_tsv(selected_genes_bulk,'Selected_gene_bulk.tsv')
write_tsv(selected_genes_sc,'Selected_gene_sc.tsv')

write_tsv(top20_gene_bulk, "Top20_marker_bulk.tsv")
write_tsv(top20_genes_sc, "Top20_marker_sc.tsv")


nnls(as.matrix(selected_genes_sc), as.vector(selected_genes_bulk[,4]))
nnls(as.matrix(top20_genes_sc), as.vector(top20_gene_bulk[,4]))

dim(selected_genes_sc)
dim(selected_genes_bulk)
head(selected_genes_bulk)
head(selected_genes_sc)
sum_to_one <- c(rep(1,19))
test <- rbind(selected_genes_sc,sum_to_one)
test[2863:2865,]
sum_to_one <- c(as.factor('constraint'),rep(1,90))
test_bulk <- rbind(selected_genes_bulk,sum_to_one)
test_bulk[2863:2865,]

nnls(as.matrix(test), as.vector(test_bulk[,4]))
