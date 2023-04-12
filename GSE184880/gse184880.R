if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
BiocManager::install("SingleR")
BiocManager::install("scRNAseq")

library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(muscat)
library(reshape2)
library(DEsingle)

setwd("~/gse184880")
a1 <- Read10X(data.dir = "~/gse184880/ca1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca1.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca2.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca3.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca4")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca4.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca5")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca5.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca6")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca6.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca7")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "ca7.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "nt1.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "nt2.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "nt3.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt4")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "nt4.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt5")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "nt5.rds")

ca1<-readRDS(file="ca1.rds")
ca2<-readRDS(file="ca2.rds")
ca3<-readRDS(file="ca3.rds")
ca4<-readRDS(file="ca4.rds")
ca5<-readRDS(file="ca5.rds")
ca6<-readRDS(file="ca6.rds")
ca7<-readRDS(file="ca7.rds")
nt1<-readRDS(file="nt1.rds")
nt2<-readRDS(file="nt2.rds")
nt3<-readRDS(file="nt3.rds")
nt4<-readRDS(file="nt4.rds")
nt5<-readRDS(file="nt5.rds")

ca1<-RenameCells(ca1,add.cell.id="ca1",for.merge=T)
ca1@meta.data$tech<-"cancer"
ca1@meta.data$celltype<-"ca1"

ca2<-RenameCells(ca2,add.cell.id="ca2",for.merge=T)
ca2@meta.data$tech<-"cancer"
ca2@meta.data$celltype<-"ca2"

ca3<-RenameCells(ca3,add.cell.id="ca3",for.merge=T)
ca3@meta.data$tech<-"cancer"
ca3@meta.data$celltype<-"ca3"

ca4<-RenameCells(ca4,add.cell.id="ca4",for.merge=T)
ca4@meta.data$tech<-"cancer"
ca4@meta.data$celltype<-"ca4"

ca5<-RenameCells(ca5,add.cell.id="ca5",for.merge=T)
ca5@meta.data$tech<-"cancer"
ca5@meta.data$celltype<-"ca5"

ca6<-RenameCells(ca6,add.cell.id="ca6",for.merge=T)
ca6@meta.data$tech<-"cancer"
ca6@meta.data$celltype<-"ca6"

ca7<-RenameCells(ca7,add.cell.id="ca7",for.merge=T)
ca7@meta.data$tech<-"cancer"
ca7@meta.data$celltype<-"ca7"

nt1<-RenameCells(nt1,add.cell.id="nt1",for.merge=T)
nt1@meta.data$tech<-"normal"
nt1@meta.data$celltype<-"nt1"

nt2<-RenameCells(nt2,add.cell.id="nt2",for.merge=T)
nt2@meta.data$tech<-"normal"
nt2@meta.data$celltype<-"nt2"

nt3<-RenameCells(nt3,add.cell.id="nt3",for.merge=T)
nt3@meta.data$tech<-"normal"
nt3@meta.data$celltype<-"nt3"

nt4<-RenameCells(nt4,add.cell.id="nt4",for.merge=T)
nt4@meta.data$tech<-"normal"
nt4@meta.data$celltype<-"nt4"

nt5<-RenameCells(nt5,add.cell.id="nt5",for.merge=T)
nt5@meta.data$tech<-"normal"
nt5@meta.data$celltype<-"nt5"

p1<-merge(ca1,nt1)
p2<-merge(ca2,nt2)
p3<-merge(ca3,nt3)
p4<-merge(ca4,nt4)
p5<-merge(ca5,nt5)
ca67<-merge(ca6,ca7)
p12<-merge(p1,p2)
p34<-merge(p3,p4)
p56<-merge(p5,ca67)
p1234<-merge(p12,p34)
p123456<-merge(p1234,p56)

saveRDS(p123456, file="p123456_before_integrate.rds")
hms<-p123456

#before integrate
#work on cancer_VS_normal HGSOC
setwd("~/gse184880/cancer_vs_normal")
hms<-readRDS("p123456_before_integrate.rds")
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
reference.list <- pancreas.list[c("ca1","ca2","ca3","ca4","ca5","ca6","ca7",
                                  "nt1","nt2","nt3","nt4","nt5")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:5)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1
#find how many 15cluster
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:10)
hms_cluster <- FindClusters( hms_neighbor, resolution = 0.5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:5)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test.rds")
cluster0.markers <- FindMarkers(hms_cluster, ident.1=0, min.pcr=0.25)
head(cluster0.markers, n=10)
cluster1.markers <- FindMarkers(hms_cluster, ident.1=1, min.pcr=0.25)
head(cluster1.markers, n=10)
cluster2.markers <- FindMarkers(hms_cluster, ident.1=2, min.pcr=0.25)
head(cluster2.markers, n=10)
cluster3.markers <- FindMarkers(hms_cluster, ident.1=3, min.pcr=0.25)
head(cluster3.markers, n=10)
cluster4.markers <- FindMarkers(hms_cluster, ident.1=4, min.pcr=0.25)
head(cluster4.markers, n=10)
cluster5.markers <- FindMarkers(hms_cluster, ident.1=5, min.pcr=0.25)
head(cluster5.markers, n=10)
cluster6.markers <- FindMarkers(hms_cluster, ident.1=6, min.pcr=0.25)
head(cluster6.markers, n=10)
cluster7.markers <- FindMarkers(hms_cluster, ident.1=7, min.pcr=0.25)
head(cluster7.markers, n=10)
cluster8.markers <- FindMarkers(hms_cluster, ident.1=8, min.pcr=0.25)
head(cluster8.markers, n=10)
cluster9.markers <- FindMarkers(hms_cluster, ident.1=9, min.pcr=0.25)
head(cluster9.markers, n=10)
cluster10.markers <- FindMarkers(hms_cluster, ident.1=10, min.pcr=0.25)
head(cluster10.markers, n=10)
cluster11.markers <- FindMarkers(hms_cluster, ident.1=11, min.pcr=0.25)
head(cluster11.markers, n=10)

new.cluster.ids <- c("NKT", "B","B", "Hepatocyte","Treg",
                     "Plasma", "NKT","Dendritic",
                     "Progenitor", "Treg","Progenitor","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Hepatocyte<-subset(hms_cluster_id, idents=c('Hepatocyte'))
DimPlot(Hepatocyte, reduction = "umap")
saveRDS(Hepatocyte, file="Hepatocyte.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

Plasma<-subset(hms_cluster_id, idents=c('Plasma'))
DimPlot(Plasma, reduction = "umap")
saveRDS(Plasma, file="Plasma.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

#input each cluster
NKT<-readRDS("NKT.rds")
B<-readRDS("B.rds")
Hepatocyte<-readRDS("Hepatocyte.rds")
Treg<-readRDS("Treg.rds")
Plasma<-readRDS("Plasma.rds")
Dendritic<-readRDS("Dendritic.rds")
Progenitor<-readRDS("Progenitor.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#deg in NKT
NKT<-readRDS("NKT.rds")
a<-NKT@meta.data
write.table(a,"a")
class(NKT)
NKT.sec<-as.SingleCellExperiment(NKT)
group<-factor(c(rep(1,2165),rep(2,363)))

rds<-readRDS('NKT.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_NKT")
write.table(results.classified,"results.classified_NKT")

#deg in Treg
Treg<-readRDS("Treg.rds")
b<-Treg@meta.data
write.table(b,"b")
class(Treg)
Treg.sec<-as.SingleCellExperiment(Treg)
group<-factor(c(rep(1,760),rep(2,17)))

rds<-readRDS('Treg.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Treg")
write.table(results.classified,"results.classified_Treg")

# early and late stage HGSOC1 IIIB, HGSOC2 IIB, HGSOC3 IC2, HGSOC4 IC2, HGSOC5 IIB, HGSOC6 IIIC, HGSOC7 IC2.
#Here, HGSOC3,4,7(IC2) as early, HGSOC1,2,5,6(IIB, IIIB and IIIC) as late. We want P value. late_VS_early
setwd("~/gse184880/late_vs_early")
ca1<-readRDS(file="ca1.rds")
ca2<-readRDS(file="ca2.rds")
ca3<-readRDS(file="ca3.rds")
ca4<-readRDS(file="ca4.rds")
ca5<-readRDS(file="ca5.rds")
ca6<-readRDS(file="ca6.rds")
ca7<-readRDS(file="ca7.rds")

ca1<-RenameCells(ca1,add.cell.id="ca1",for.merge=T)
ca1@meta.data$tech<-"late"
ca1@meta.data$celltype<-"ca1"

ca2<-RenameCells(ca2,add.cell.id="ca2",for.merge=T)
ca2@meta.data$tech<-"late"
ca2@meta.data$celltype<-"ca2"

ca3<-RenameCells(ca3,add.cell.id="ca3",for.merge=T)
ca3@meta.data$tech<-"early"
ca3@meta.data$celltype<-"ca3"

ca4<-RenameCells(ca4,add.cell.id="ca4",for.merge=T)
ca4@meta.data$tech<-"early"
ca4@meta.data$celltype<-"ca4"

ca5<-RenameCells(ca5,add.cell.id="ca5",for.merge=T)
ca5@meta.data$tech<-"late"
ca5@meta.data$celltype<-"ca5"

ca6<-RenameCells(ca6,add.cell.id="ca6",for.merge=T)
ca6@meta.data$tech<-"late"
ca6@meta.data$celltype<-"ca6"

ca7<-RenameCells(ca7,add.cell.id="ca7",for.merge=T)
ca7@meta.data$tech<-"early"
ca7@meta.data$celltype<-"ca7"

ca12<-merge(ca1,ca2)
ca34<-merge(ca3,ca4)
ca56<-merge(ca5,ca6)
ca1256<-merge(ca12,ca56)
ca347<-merge(ca34,ca7)
ca1256347<-merge(ca1256,ca347)

saveRDS(ca1256347, file="ca1256347_before_integrate.rds")
hms<-ca1256347

#before integrate
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
reference.list <- pancreas.list[c("ca1","ca2","ca5","ca6","ca3","ca4","ca7"
                                  )]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:5)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1
#find how many 15cluster
ElbowPlot(hms_individual_integrated)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:10)
hms_cluster <- FindClusters( hms_neighbor, resolution = 0.5)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:5)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test.rds")

hms_cluster<-readRDS(file="hms_cluster_test.rds")
cluster0.markers <- FindMarkers(hms_cluster, ident.1=0, min.pcr=0.25)
head(cluster0.markers, n=10)
cluster1.markers <- FindMarkers(hms_cluster, ident.1=1, min.pcr=0.25)
head(cluster1.markers, n=10)
cluster2.markers <- FindMarkers(hms_cluster, ident.1=2, min.pcr=0.25)
head(cluster2.markers, n=10)
cluster3.markers <- FindMarkers(hms_cluster, ident.1=3, min.pcr=0.25)
head(cluster3.markers, n=10)
cluster4.markers <- FindMarkers(hms_cluster, ident.1=4, min.pcr=0.25)
head(cluster4.markers, n=10)
cluster5.markers <- FindMarkers(hms_cluster, ident.1=5, min.pcr=0.25)
head(cluster5.markers, n=10)
cluster6.markers <- FindMarkers(hms_cluster, ident.1=6, min.pcr=0.25)
head(cluster6.markers, n=10)
cluster7.markers <- FindMarkers(hms_cluster, ident.1=7, min.pcr=0.25)
head(cluster7.markers, n=10)
cluster8.markers <- FindMarkers(hms_cluster, ident.1=8, min.pcr=0.25)
head(cluster8.markers, n=10)
cluster9.markers <- FindMarkers(hms_cluster, ident.1=9, min.pcr=0.25)
head(cluster9.markers, n=10)
cluster10.markers <- FindMarkers(hms_cluster, ident.1=10, min.pcr=0.25)
head(cluster10.markers, n=10)
cluster11.markers <- FindMarkers(hms_cluster, ident.1=11, min.pcr=0.25)
head(cluster11.markers, n=10)

new.cluster.ids <- c("NKT", "B","B", "Hepatocyte","Treg",
                     "Plasma", "NKT","Dendritic",
                     "Progenitor", "Treg","Progenitor","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Hepatocyte<-subset(hms_cluster_id, idents=c('Hepatocyte'))
DimPlot(Hepatocyte, reduction = "umap")
saveRDS(Hepatocyte, file="Hepatocyte.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

Plasma<-subset(hms_cluster_id, idents=c('Plasma'))
DimPlot(Plasma, reduction = "umap")
saveRDS(Plasma, file="Plasma.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

#input each cluster
NKT<-readRDS("NKT.rds")
B<-readRDS("B.rds")
Hepatocyte<-readRDS("Hepatocyte.rds")
Treg<-readRDS("Treg.rds")
Plasma<-readRDS("Plasma.rds")
Dendritic<-readRDS("Dendritic.rds")
Progenitor<-readRDS("Progenitor.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#deg in NKT
NKT<-readRDS("NKT.rds")
a<-NKT@meta.data
write.table(a,"a")
class(NKT)
NKT.sec<-as.SingleCellExperiment(NKT)
group<-factor(c(rep(1,2165),rep(2,363)))

rds<-readRDS('NKT.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_NKT")
write.table(results.classified,"results.classified_NKT")

#deg in Treg
Treg<-readRDS("Treg.rds")
b<-Treg@meta.data
write.table(b,"b")
class(Treg)
Treg.sec<-as.SingleCellExperiment(Treg)
group<-factor(c(rep(1,760),rep(2,17)))

rds<-readRDS('Treg.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Treg")
write.table(results.classified,"results.classified_Treg")

#late and early GSE184880 HGSOC

setwd("~/gse184880/late_vs_early")
hms_cluster<-readRDS(file="hms_cluster_test.rds")

new.cluster.ids <- c("CD8", "Treg","B", "NKT","B","Treg",
                     "Cancer cell", "Macrophage","Progenitor",
                     "Astrocyte", "Dendritic","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

CD8<-subset(hms_cluster_id, idents=c('CD8'))
DimPlot(CD8, reduction = "umap")
saveRDS(CD8, file="CD8.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

NKT<-subset(hms_cluster_id, idents=c('NKT'))
DimPlot(NKT, reduction = "umap")
saveRDS(NKT, file="NKT.rds")

Cancer_cell<-subset(hms_cluster_id, idents=c('Cancer cell'))
DimPlot(Cancer_cell, reduction = "umap")
saveRDS(Cancer_cell, file="Cancer_cell.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Progenitor<-subset(hms_cluster_id, idents=c('Progenitor'))
DimPlot(Progenitor, reduction = "umap")
saveRDS(Progenitor, file="Progenitor.rds")

Astrocyte<-subset(hms_cluster_id, idents=c('Astrocyte'))
DimPlot(Astrocyte, reduction = "umap")
saveRDS(Astrocyte, file="Astrocyte.rds")

Dendritic<-subset(hms_cluster_id, idents=c('Dendritic'))
DimPlot(Dendritic, reduction = "umap")
saveRDS(Dendritic, file="Dendritic.rds")



#deg in NKT
NKT<-readRDS("NKT.rds")
a<-NKT@meta.data
write.table(a,"a")
class(NKT)
NKT.sec<-as.SingleCellExperiment(NKT)
group<-factor(c(rep(1,396),rep(2,306)))

rds<-readRDS('NKT.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_NKT")
write.table(results.classified,"results.classified_NKT")

#deg in Treg
Treg<-readRDS("Treg.rds")
b<-Treg@meta.data
write.table(b,"b")
class(Treg)
Treg.sec<-as.SingleCellExperiment(Treg)
group<-factor(c(rep(1,723),rep(2,670)))

rds<-readRDS('Treg.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results_Treg,"results_Treg")
write.table(results.classified,"results.classified_Treg")

#deg in CD8
CD8<-readRDS("CD8.rds")
a<-CD8@meta.data
write.table(a,"a")
class(CD8)
CD8.sec<-as.SingleCellExperiment(NKT)
group<-factor(c(rep(1,736),rep(2,893)))

rds<-readRDS('CD8.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results_CD8,"results_CD8")
write.table(results.classified,"results.classified_CD8")
