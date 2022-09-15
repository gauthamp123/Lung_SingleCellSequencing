# install packages Seurat, dplyr and patchwork
install.packages("BiocManager")
BiocManager::install("multtest")

install.packages('Seurat')
library(Seurat)
library(patchwork)
library(dplyr)

# load data
# pbmc_ser_data <- Read10X(data.dir = ".../curr_dir/filtered_gene_bc_matrices/hg19")
pbmc_ser_data <- Read10X(data.dir = "/Users/priyu/OneDrive/Documents/R/CSE185/filtered_gene_bc_matrices/hg19")

pbmc_ser <- CreateSeuratObject(counts = pbmc_ser_data, project = "pbmc_ser3k", min.cells = 3, min.features = 200)

pbmc_ser[["percent.mt"]] <- PercentageFeatureSet(pbmc_ser, pattern = "^MT-")

pbmc_ser <- subset(pbmc_ser, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc_ser <- NormalizeData(pbmc_ser, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_ser <- NormalizeData(pbmc_ser)

pbmc_ser <- FindVariableFeatures(pbmc_ser, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc_ser), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc_ser)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


all.genes <- rownames(pbmc_ser)
pbmc_ser <- ScaleData(pbmc_ser, features = all.genes)

pbmc_ser <- RunPCA(pbmc_ser, features = VariableFeatures(object = pbmc_ser))

# Find clusters
pbmc_ser <- FindNeighbors(pbmc_ser, dims = 1:10)
pbmc_ser <- FindClusters(pbmc_ser, resolution = 0.5)

# Create the pca graph
pbmc_ser <- RunPCA(pbmc_ser, features = VariableFeatures(object = pbmc_ser))
DimPlot(pbmc_ser, reduction = "pca")

# Create tSNE
pbmc_ser <- RunTSNE(pbmc_ser, dims = 1:10)
DimPlot(pbmc_ser, reduction = "tsne")

FeaturePlot(pbmc_ser, features = c("CD3E", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
