setwd("~/gse184880/cancer_vs_normal")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="T", "1"="Fibroblast","2"="Epithelial","3"= "T",
                                                 "4"="B", "5"="Macrophage","6"="Macrophage",
                                                 "7"="Endothelial", "8"="B","9"="B")

colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

plotCytoTRACE(result01, phenotype = phe)

setwd("~/gse184880/cancer_vs_normal")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="t.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="CD8_T",  "1"="CD8_T","2"="CD4_T",  "3"="CD8_T","4"="CD4_T","5"="NK",
                                                 "6"="Treg", "7"="Treg", "8"="CD8_ pre_exhausted_T", "9"="CD8_cytotoxicity_T",
                                                 "10"="CD8_cytotoxicity_T","11"="CD8_T")

colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

plotCytoTRACE(result01, phenotype = phe)

setwd("~/gse184880/late_vs_early")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="T", "1"="B","2"="Epithelial","3"= "T","4"="T","5"="Progenitor","6"="T",
                                                 "7"="Macrophage","8"="Macrophage","9"="T","10"="NKT","11"="T","12"="Epithelial",
                                                 "13"="Macrophage","14"="Endothelial","15"="Endothelial","16"="B","17"="Dendritic")

colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

plotCytoTRACE(result01, phenotype = phe)

setwd("~/gse184880/late_vs_early")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="t.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="CD8_cytotoxic_T",  "1"="NK","2"="Treg",  "3"="CD4_effector_memory_T",
                                                 "4"="CD4_transitional_memory_T","5"="CD8_cytotoxic_T","6"="CD8_tissue_resident_memory_T",
                                                 "7"="CD8_terminally_exhausted_T", "8"="CD8_pre_exhausted_T", "9"="B","10"="Treg",
                                                 "11"="CD8_pre_exhausted_T","12"="B","13"="Epithelial","14"="Conventional_natural_killer")

colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)

plotCytoTRACE(result01, phenotype = phe)
