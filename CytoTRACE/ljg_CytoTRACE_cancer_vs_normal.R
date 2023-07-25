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
a1 <- Read10X(data.dir = "~/gse184880/ca1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
#ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca1.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca2.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca3.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca4")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca4.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca5")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca5.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca6")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca6.rds")

a1 <- Read10X(data.dir = "~/gse184880/ca7")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "ca7.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt1")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "nt1.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt2")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "nt2.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt3")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "nt3.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt4")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
saveRDS(pbmc, file = "nt4.rds")

a1 <- Read10X(data.dir = "~/gse184880/nt5")
pbmc <- CreateSeuratObject(counts = a1, project = "a1", min.cells =3, min.features=200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#查看nFeature_RNA
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#把原代碼中"nFeature_RNA < 2500 &"刪去 ,下面是源代碼
#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:20)
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

ca12<-merge(ca1,ca2)
ca34<-merge(ca3,ca4)
ca56<-merge(ca5,ca6)
ca1234<-merge(ca12,ca34)
ca567<-merge(ca56,ca7)
ca1234567<-merge(ca1234,ca567)

nt12<-merge(nt1,nt2)
nt34<-merge(nt3,nt4)
nt345<-merge(nt34,nt5)
nt12345<-merge(nt12,nt345)

ca_nt<-merge(ca1234567,nt12345)
saveRDS(ca_nt, file="ca_nt_before_integrate.rds")

#before integrate
#work on cancer_VS_normal HGSOC
setwd("~/gse184880")
hms<-readRDS("ca_nt_before_integrate.rds")
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
hms_cluster <- FindClusters(hms_neighbor, resolution = 0.2)#第1個聚類只有0個基因,調整resolution的值,往下調
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

new.cluster.ids <- c("T", "Fibroblast","Epithelial", "T",
                     "B", "Macrophage","Macrophage",
                     "Endothelial", "B","B") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_test.rds")

hms_cluster_id<-readRDS(file="hms_cluster_id_test.rds")
T<-subset(hms_cluster_id, idents=c('T'))
DimPlot(T, reduction = "umap")
saveRDS(T, file="T.rds")

B<-subset(hms_cluster_id, idents=c('B'))
DimPlot(B, reduction = "umap")
saveRDS(B, file="B.rds")

Fibroblast<-subset(hms_cluster_id, idents=c('Fibroblast'))
DimPlot(Fibroblast, reduction = "umap")
saveRDS(Fibroblast, file="Fibroblast.rds")

Epithelial<-subset(hms_cluster_id, idents=c('Epithelial'))
DimPlot(Epithelial, reduction = "umap")
saveRDS(Epithelial, file="Epithelial.rds")

Macrophage<-subset(hms_cluster_id, idents=c('Macrophage'))
DimPlot(Macrophage, reduction = "umap")
saveRDS(Macrophage, file="Macrophage.rds")

Endothelial<-subset(hms_cluster_id, idents=c('Endothelial'))
DimPlot(Endothelial, reduction = "umap")
saveRDS(Endothelial, file="Endothelial.rds")

#input each cluster
T<-readRDS("T.rds")
B<-readRDS("B.rds")
Fibroblast<-readRDS("Fibroblast.rds")
Epithelial<-readRDS("Epithelial.rds")
Macrophage<-readRDS("Macrophage.rds")
Endothelial<-readRDS("Endothelial.rds")
hms_cluster_id<-readRDS("hms_cluster_id_test.rds")

#regroup T cell
T<-readRDS(file="T.rds")
DimPlot(T, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dim(T@meta.data)
#T <- NormalizeData(T, normalization.method = "LogNormalize", scale.factor = 10000)
T <- FindVariableFeatures(T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(T)
#T- ScaleData(T, rownames(T))#跑不出來
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

setwd("~/gse184880")
hms_cluster<-readRDS(file="T_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap")

new.cluster.ids <- c("CD8_T","CD8_T","CD4_T","CD8_T","CD4_T","Natural_killer",
                     "Treg","Treg","CD8_pre_exhausted_T","CD8_cytotoxicity_T",
                     "CD8_cytotoxicity_T", "CD8_T") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
names(hms_cluster_id@meta.data)

DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "T_cluster_id_test.rds")

setwd("~/gse184880")
hms_cluster_id<-readRDS(file="T_cluster_id_test.rds")

#each type of cells
CD4_T<-subset(hms_cluster_id, idents=c('CD4_T'))
DimPlot(CD4_T, reduction = "umap")
saveRDS(CD4_T, file="CD4_T.rds")

Treg<-subset(hms_cluster_id, idents=c('Treg'))
DimPlot(Treg, reduction = "umap")
saveRDS(Treg, file="Treg.rds")

CD8_T<-subset(hms_cluster_id, idents=c('CD8_T'))
DimPlot(CD8_T, reduction = "umap")
saveRDS(CD8_T, file="CD8_T.rds")

Natural_killer<-subset(hms_cluster_id, idents=c('Natural_killer'))
DimPlot(Natural_killer, reduction = "umap")
saveRDS(Natural_killer, file="Natural_killer.rds")

CD8_pre_exhausted_T<-subset(hms_cluster_id, idents=c('CD8_pre_exhausted_T'))
DimPlot(CD8_pre_exhausted_T, reduction = "umap")
saveRDS(CD8_pre_exhausted_T, file="CD8_pre_exhausted_T.rds")

CD8_cytotoxicity_T<-subset(hms_cluster_id, idents=c('CD8_cytotoxicity_T'))
DimPlot(CD8_cytotoxicity_T, reduction = "umap")
saveRDS(CD8_cytotoxicity_T, file="CD8_cytotoxicity_T.rds")

#新增代码：去除非T细胞
T_new <- subset(hms_cluster_id, idents = c("CD8_T","CD4_T","Treg",
                                           "CD8_pre_exhausted_T","CD8_cytotoxicity_T"))
saveRDS(T_new, file = "t.rds")
DimPlot(T_new, reduction = "umap", label = FALSE, pt.size = 0.5) 
DimPlot(T_new, reduction = "umap", label = TRUE, pt.size =0.5) 

setwd("~/gse184880")
hms_cluster<-readRDS(file="t.rds")
DimPlot(hms_cluster, reduction = "umap")

VlnPlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R"))
VlnPlot(hms_cluster, features = c("CD4","FOXP3","CTLA4","TRDC","TGFB1","TRAV1-2"))
VlnPlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2", "IGHG1","MZB1"))

RidgePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R"))
RidgePlot(hms_cluster, features = c("CD4","FOXP3","CTLA4","TRDC","TGFB1","TRAV1-2"))
RidgePlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2", "IGHG1","MZB1"))

FeaturePlot(hms_cluster, features = c("CD8A", "CD8B","GZMK","CD27","CCR7","IL7R"))
FeaturePlot(hms_cluster, features = c("CD4","FOXP3","CTLA4","TRDC","TGFB1","TRAV1-2"))
FeaturePlot(hms_cluster, features = c("LAG3","PDCD1","TGFBR2", "IGHG1","MZB1"))

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
group<-factor(c(rep(1,730),rep(2,14)))

rds<-readRDS('Treg.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_Treg")
write.table(results.classified,"results.classified_Treg")

#deg in CD8_T
CD8_T<-readRDS("CD8_T.rds")
b<-CD8_T@meta.data
write.table(b,"b")
class(CD8_T)
CD8_T.sec<-as.SingleCellExperiment(CD8_T)
group<-factor(c(rep(1,1629),rep(2,126)))

rds<-readRDS('CD8_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_T")
write.table(results.classified,"results.classified_CD8_T")

#deg in CD4_T
CD4_T<-readRDS("CD4_T.rds")
c<-CD4_T@meta.data
write.table(c,"c")
class(CD4_T)
CD4_T.sec<-as.SingleCellExperiment(CD4_T)
group<-factor(c(rep(1,879),rep(2,85)))

rds<-readRDS('CD4_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD4_T")
write.table(results.classified,"results.classified_CD4_T")

#deg in CD8_pre_exhausted_T
CD8_pre_exhausted_T<-readRDS("CD8_pre_exhausted_T.rds")
d<-CD8_pre_exhausted_T@meta.data
write.table(d,"d")
class(CD8_pre_exhausted_T)
CD8_pre_exhausted_T.sec<-as.SingleCellExperiment(CD8_pre_exhausted_T)
group<-factor(c(rep(1,204),rep(2,1)))

rds<-readRDS('CD8_pre_exhausted_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_pre_exhausted_T")
write.table(results.classified,"results.classified_CD8_pre_exhausted_T")

#deg in CD8_cytotoxicity_T
CD8_cytotoxicity_T<-readRDS("CD8_cytotoxicity_T.rds")
e<-CD8_cytotoxicity_T@meta.data
write.table(e,"e")
class(CD8_cytotoxicity_T)
CD8_cytotoxicity_T.sec<-as.SingleCellExperiment(CD8_cytotoxicity_T)
group<-factor(c(rep(1,292),rep(2,58)))

rds<-readRDS('CD8_cytotoxicity_T.rds')
counts<-as.matrix(rds@assays$RNA@counts)
results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results,"results_CD8_cytotoxicity_T")
write.table(results.classified,"results.classified_CD8_cytotoxicity_T")

#volcano_plot
setwd("~/gse184880/cancer_vs_normal/volcano_data")
res <- read.csv("cd4t_volcano.txt", header=TRUE,sep="\t")
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

res <- read.csv("cd8_cytotoxicity_t.txt", header=TRUE,sep="\t")
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

res <- read.csv("cd8t_volcano.txt", header=TRUE,sep="\t")
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

res <- read.csv("treg_volcano.txt", header=TRUE,sep="\t")
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