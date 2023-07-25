setwd("~/gse184880/cancer_vs_normal")
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/
#library(monocle)#can not use this package.error will come out.
#cancer_vs_normal all cells
seurat <- readRDS(file="hms_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis

#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4

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

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
 "0"="T", "1"="Fibroblast","2"="Epithelial","3"= "T",
 "4"="B", "5"="Macrophage","6"="Macrophage",
 "7"="Endothelial", "8"="B","9"="B" )
colnames(colData(cds))

cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2
p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3
p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 6)
p5

##细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)

##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
#long time
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)#take long time

write.csv(Track_genes,"Trajectory_genes_all.csv",row.names = TRUE)

#20 genes    morans_I>0.6 compared with 434_immune_gene
Track<-c("FOXP3","IL32","CD8A","CD3E","FCRL5","CD3G","TNFRSF17","GZMA","PTPRC","CD79A",
         "CD52","CD2","CCL5","CD3D")

plot_genes_in_pseudotime(cds[Track,], color_cells_by="tech", 
                         min_expr=0.5, ncol = 4)

plot_genes_in_pseudotime(cds[Track,], color_cells_by="assigned_cell_type", 
                         min_expr=0.5, ncol = 4)


#FeaturePlot图
plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)


#cancer_vs_normal only T cells
seurat <- readRDS(file="t.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis

#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4

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

colData(cds)$assigned_cell_type <- as.character(partitions(cds))

colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
 "0"="CD8_T",  "1"="CD8_T","2"="CD4_T",  "3"="CD8_T","4"="CD4_T","5"="NK",
 "6"="Treg", "7"="Treg", "8"="CD8_ pre_exhausted_T", "9"="CD8_cytotoxicity_T",
 "10"="CD8_cytotoxicity_T","11"="CD8_T" )
colnames(colData(cds))

cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2
p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3
p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 6)
p5

##细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)

##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
#long time
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)#take long time

write.csv(Track_genes,"Trajectory_genes.csv",row.names = TRUE)

#20 genes    morans_I>0.5 compared with 434_immune_gene
Track<-c("CCL5","FOXP3","GZMA","CD8A","TNFRSF4","TNFRSF18","IL7R","GZMH",
"CD8B","GZMB","IL2RA")

plot_genes_in_pseudotime(cds[Track,], color_cells_by="tech", 
                         min_expr=0.5, ncol = 4)

plot_genes_in_pseudotime(cds[Track,], color_cells_by="assigned_cell_type", 
                         min_expr=0.5, ncol = 4)
#FeaturePlot图
plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)