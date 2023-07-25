if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
BiocManager::install("SingleR")
BiocManager::install("scRNAseq")

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(dplyr)
library(Matrix)
library(reshape2)
library(celldex)
library(clustree)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(muscat)
library(DEsingle)

setwd("~/gse184880")
ca1<-readRDS(file="ca1.rds")
ca2<-readRDS(file="ca2.rds")
ca3<-readRDS(file="ca3.rds")
ca4<-readRDS(file="ca4.rds")
ca5<-readRDS(file="ca5.rds")
ca6<-readRDS(file="ca6.rds")
ca7<-readRDS(file="ca7.rds")

ca1<-RenameCells(ca1,add.cell.id="ca1",for.merge=T)
ca1@meta.data$tech<-"IIIB"
ca1@meta.data$celltype<-"ca1"

ca2<-RenameCells(ca2,add.cell.id="ca2",for.merge=T)
ca2@meta.data$tech<-"IIB"
ca2@meta.data$celltype<-"ca2"

ca3<-RenameCells(ca3,add.cell.id="ca3",for.merge=T)
ca3@meta.data$tech<-"IC2"
ca3@meta.data$celltype<-"ca3"

ca4<-RenameCells(ca4,add.cell.id="ca4",for.merge=T)
ca4@meta.data$tech<-"IC2"
ca4@meta.data$celltype<-"ca4"

ca5<-RenameCells(ca5,add.cell.id="ca5",for.merge=T)
ca5@meta.data$tech<-"IIB"
ca5@meta.data$celltype<-"ca5"

ca6<-RenameCells(ca6,add.cell.id="ca6",for.merge=T)
ca6@meta.data$tech<-"IIIC"
ca6@meta.data$celltype<-"ca6"

ca7<-RenameCells(ca7,add.cell.id="ca7",for.merge=T)
ca7@meta.data$tech<-"IC2"
ca7@meta.data$celltype<-"ca7"

ca61<-merge(ca6,ca1)
ca25<-merge(ca2,ca5)
ca6125<-merge(ca61,ca25)
ca34<-merge(ca3,ca4)
ca347<-merge(ca34,ca7)
ca<-merge(ca6125,ca347)
saveRDS(ca, file="ca_before_integrate.rds")

#before integrate,work on late_VS_eary HGSOC 
## early and late stage HGSOC1 IIIB, HGSOC2 IIB, HGSOC3 IC2, HGSOC4 IC2, HGSOC5 IIB, HGSOC6 IIIC, HGSOC7 IC2.
#Here, HGSOC2,3,4,5,7(IC2,IIB) as early, HGSOC6,1(IIIB and IIIC) as late.  late_VS_early
setwd("~/gse184880/late_vs_early")
hms<-readRDS("ca_before_integrate.rds")
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) 
plot_grid(p1,p2)
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}
reference.list <- pancreas.list[c("ca6","ca1","ca2","ca5","ca3","ca4","ca7")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:20)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1
#find how many 15cluster
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:20)
#obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
#clustree(obj)
hms_cluster <- FindClusters(hms_neighbor, resolution = 0.8)#第1個聚類只有0個基因,調整resolution的值,往下調
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:20)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.25, #如果設定是0.25,第1個聚類只有0個基因,無法write.csv
                                logfc.threshold = 0.25
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(top20_table,"top20_marker_genes.csv",row.names=FALSE)

#细胞及细胞中基因与RNA数量
slotNames(hms)
#assay
hms@assays
dim(hms@meta.data)
View(hms@meta.data)

hms_cluster<-readRDS(file = "hms_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")
new.cluster.ids <- c("T", "B","Epithelial", "T","T","Progenitor","T","Macrophage",
                     "Macrophage","T","NKT","T","Epithelial","Macrophage",
                     "Endothelial","Endothelial", "B","Dendritic") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

setwd("~/gse184880/late_vs_early")
hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
T<-subset(hms_cluster_id, idents=c('T'))
DimPlot(T, reduction = "umap")
saveRDS(T, file="T.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")


#input each cluster
T<-readRDS("T.rds")
B<-readRDS("B.rds")
Epithelial<-readRDS("Epithelial.rds")
Progenitor<-readRDS("Progenitor.rds")
Macrophage<-readRDS("Macrophage.rds")
NKT<-readRDS("NKT.rds")
Endothelial<-readRDS("Endothelial.rds")
Dendritic<-readRDS("Dendritic.rds")
hms_cluster_id<-readRDS("hms_cluster_id_test.rds")


#regroup T cell
T_new <- subset(hms_cluster_id, idents = c("T","NKT"))
saveRDS(T_new, file = "t.rds")
DimPlot(T_new, reduction = "umap", label = FALSE, pt.size = 0.5) 
DimPlot(T_new, reduction = "umap", label = TRUE, pt.size =0.5) 

T<-readRDS(file="t.rds")
VlnPlot(T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(T@meta.data)
#T <- NormalizeData(T, normalization.method = "LogNormalize", scale.factor = 10000)
T <- FindVariableFeatures(T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(T)
T<- ScaleData(T, rownames(T))#跑不出來
T <- RunPCA(T, features = VariableFeatures(object = T))
#Determine the ‘dimensionality’ of the dataset
T <- JackStraw(T, num.replicate = 100,dims = 40)
T <- ScoreJackStraw(T, dims = 1:40)
JackStrawPlot(T, dims = 1:40)
ElbowPlot(T,ndims = 40)
#確定下面的dims
hms_neighbor<- FindNeighbors(T, dims = 1:20)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
clustree(obj)
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
hms_cluster <- FindClusters(hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:10)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "T_cluster_test.rds")
hms_cluster<-readRDS(file="T_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")
scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.1, #如果設定是0.25,第11個聚類只有14個基因,無法write.csv
                                logfc.threshold = 0.25
)
write.table(scRNA.markers,file="T_cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="T_top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(top20_table,"T_top20_marker_genes.csv",row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(T)
#assay
T@assays
dim(T@meta.data)
View(T@meta.data)

setwd("~/gse184880/late_vs_early")
hms_cluster<-readRDS(file="T_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c("CD8_cytotoxic_T","Natural_killer","Treg","CD4_effector_memory_T",
                     "CD4_transitional_memory_T","CD8_cytotoxic_T",
                     "CD8_tissue_resident_memory_T", "CD8_terminally_exhausted_T",
                     "CD8_pre_exhausted_T", "B","Treg","CD8_pre_exhausted_T","B",
                     "Epithelial", "Conventional_natural_killer") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
names(hms_cluster_id@meta.data)

DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "T_cluster_id_test.rds")

setwd("~/gse184880/late_vs_early")
hms_cluster_id<-readRDS(file="T_cluster_id_test.rds")

#each type of cells
CD8_cytotoxic_T<-subset(hms_cluster_id, idents=c('CD8_cytotoxic_T'))
DimPlot(CD8_cytotoxic_T, reduction = "umap")
saveRDS(CD8_cytotoxic_T, file="CD8_cytotoxic_T.rds")

Natural_killer<-subset(hms_cluster_id, idents=c('Natural_killer'))
DimPlot(Natural_killer, reduction = "umap")
saveRDS(Natural_killer, file="Natural_killer.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

CD4_effector_memory_T<-subset(hms_cluster_id, idents=c('CD4_effector_memory_T'))
DimPlot(CD4_effector_memory_T, reduction = "umap")
saveRDS(CD4_effector_memory_T, file="CD4_effector_memory_T.rds")

CD4_transitional_memory_T<-subset(hms_cluster_id, idents=c('CD4_transitional_memory_T'))
DimPlot(CD4_transitional_memory_T, reduction = "umap")
saveRDS(CD4_transitional_memory_T, file="CD4_transitional_memory_T.rds")

CD8_tissue_resident_memory_T<-subset(hms_cluster_id, idents=c('CD8_tissue_resident_memory_T'))
DimPlot(CD8_tissue_resident_memory_T, reduction = "umap")
saveRDS(CD8_tissue_resident_memory_T, file="CD8_tissue_resident_memory_T.rds")

CD8_terminally_exhausted_T<-subset(hms_cluster_id, idents=c('CD8_terminally_exhausted_T'))
DimPlot(CD8_terminally_exhausted_T, reduction = "umap")
saveRDS(CD8_terminally_exhausted_T, file="CD8_terminally_exhausted_T.rds")

CD8_pre_exhausted_T<-subset(hms_cluster_id, idents=c('CD8_pre_exhausted_T'))
DimPlot(CD8_pre_exhausted_T, reduction = "umap")
saveRDS(CD8_pre_exhausted_T, file="CD8_pre_exhausted_T.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

Conventional_natural_killer<-subset(hms_cluster_id, idents=c('Conventional_natural_killer'))
DimPlot(Conventional_natural_killer, reduction = "umap")
saveRDS(Conventional_natural_killer, file="Conventional_natural_killer.rds")

#新增代码：去除非T细胞
T_new <- subset(hms_cluster_id, idents = c("CD8_cytotoxic_T","Treg","CD4_effector_memory_T",
 "CD4_transitional_memory_T","CD8_tissue_resident_memory_T","CD8_terminally_exhausted_T",
                                           "CD8_pre_exhausted_T"))

saveRDS(T_new, file = "t.rds")
DimPlot(T_new, reduction = "umap", label = FALSE, pt.size = 0.5) 
DimPlot(T_new, reduction = "umap", label = TRUE, pt.size =0.5) 

setwd("~/gse184880/late_vs_early")
hms_cluster<-readRDS(file="t.rds")
DimPlot(hms_cluster, reduction = "umap")

VlnPlot(hms_cluster, features = c("CD8A", "CD8B","GZMK"))
VlnPlot(hms_cluster, features = c("CD27","CCR7","IL7R"))
VlnPlot(hms_cluster, features = c("CD4","FOXP3","CTLA4"))
VlnPlot(hms_cluster, features = c("TRDC","TGFB1","TRAV1-2"))
VlnPlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2"))
VlnPlot(hms_cluster, features = c("IGHG1","MZB1"))

RidgePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R"))
RidgePlot(hms_cluster, features = c("CD4","FOXP3","CTLA4","TRDC","TGFB1","TRAV1-2"))
RidgePlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2", "IGHG1","MZB1"))

FeaturePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R"),ncol = 3)
FeaturePlot(hms_cluster, features = c("CD4","FOXP3","CTLA4","TRDC","TGFB1","TRAV1-2"),ncol = 3)
FeaturePlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2", "IGHG1","MZB1"),ncol = 3)

#ridgeplot
RidgePlot(hms_cluster, features =c("TTN"))

#feature plot
FeaturePlot(hms_cluster, features = c("GBP5"))

#deg in Treg
Treg<-readRDS("Treg.rds")
a<-Treg@meta.data
write.table(a,"a")
class(Treg)
Treg.sec<-as.SingleCellExperiment(Treg)
group<-factor(c(rep(1,73),rep(2,609)))

rds<-readRDS('Treg.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Treg")
write.table(results.classified,"results.classified_Treg")

#deg in CD8_cytotoxic_T
CD8_cytotoxic_T<-readRDS("CD8_cytotoxic_T.rds")
b<-CD8_cytotoxic_T@meta.data
write.table(b,"b")
class(CD8_cytotoxic_T)
CD8_cytotoxic_T.sec<-as.SingleCellExperiment(CD8_cytotoxic_T)
group<-factor(c(rep(1,135),rep(2,832)))

rds<-readRDS('CD8_cytotoxic_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_cytotoxic_T.rds")
write.table(results.classified,"results.classified_CD8_cytotoxic_T")

#deg in CD8_pre_exhausted_T
CD8_pre_exhausted_T<-readRDS("CD8_pre_exhausted_T.rds")
c<-CD8_pre_exhausted_T@meta.data
write.table(c,"c")
class(CD8_pre_exhausted_T)
CD8_pre_exhausted_T.sec<-as.SingleCellExperiment(CD8_pre_exhausted_T)
group<-factor(c(rep(1,90),rep(2,366)))

rds<-readRDS('CD8_pre_exhausted_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_pre_exhausted_T")
write.table(results.classified,"results.classified_CD8_pre_exhausted_T")

#deg in CD8_tissue_resident_memory_T
CD8_tissue_resident_memory_T<-readRDS("CD8_tissue_resident_memory_T.rds")
d<-CD8_tissue_resident_memory_T@meta.data
write.table(d,"d")
class(CD8_tissue_resident_memory_T)
CD8_tissue_resident_memory_T.sec<-as.SingleCellExperiment(CD8_tissue_resident_memory_T)
group<-factor(c(rep(1,60),rep(2,286)))

rds<-readRDS('CD8_tissue_resident_memory_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_tissue_resident_memory_T")
write.table(results.classified,"results.classified_CD8_tissue_resident_memory_T")

#deg in CD8_terminally_exhausted_T
CD8_terminally_exhausted_T<-readRDS("CD8_terminally_exhausted_T.rds")
e<-CD8_terminally_exhausted_T@meta.data
write.table(e,"e")
class(CD8_terminally_exhausted_T)
CD8_terminally_exhausted_T.sec<-as.SingleCellExperiment(CD8_terminally_exhausted_T)
group<-factor(c(rep(1,5),rep(2,268)))

rds<-readRDS('CD8_terminally_exhausted_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_terminally_exhausted_T")
write.table(results.classified,"results.classified_CD8_terminally_exhausted_T")

#deg in CD4_effector_memory_T
CD4_effector_memory_T<-readRDS("CD4_effector_memory_T.rds")
f<-CD4_effector_memory_T@meta.data
write.table(f,"f")
class(CD4_effector_memory_T)
CD4_effector_memory_T.sec<-as.SingleCellExperiment(CD4_effector_memory_T)
group<-factor(c(rep(1,47),rep(2,373)))

rds<-readRDS('CD4_effector_memory_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD4_effector_memory_T")
write.table(results.classified,"results.classified_CD4_effector_memory_T")

#deg in CD4_effector_memory_T
CD4_transitional_memory_T<-readRDS("CD4_transitional_memory_T.rds")
g<-CD4_transitional_memory_T@meta.data
write.table(g,"g")
class(CD4_transitional_memory_T)
CD4_transitional_memory_T.sec<-as.SingleCellExperiment(CD4_transitional_memory_T)
group<-factor(c(rep(1,8),rep(2,411)))

rds<-readRDS('CD4_transitional_memory_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD4_transitional_memory_T")
write.table(results.classified,"results.classified_CD4_transitional_memory_T")

#volcano_plot
setwd("~/gse184880/late_vs_early/volcano_data")
res <- read.csv("CD4_effector_memory_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("CD4_transitional_memory_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-9,5),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("CD8_cytotoxic_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-7,7),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("CD8_pre_exhausted_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,11),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("CD8_terminally_exhausted_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,6),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("CD8_tissue_resident_memory_T.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,8),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

res <- read.csv("treg.txt", header=TRUE,sep="\t")
head(res)
with(res, plot(log2FoldChange, -log10(fdr), pch=20, main="Volcano plot", xlim=c(-6,8),col="grey"))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
#with(subset(res, fdr<.05 ), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="orange"))
with(subset(res, fdr<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(fdr), pch=20, col="blue"))
with(subset(res, fdr<.05 & log2FoldChange>1), points(log2FoldChange, -log10(fdr), pch=20, col="red"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

