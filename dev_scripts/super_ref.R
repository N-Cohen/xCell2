library(tidyverse)
library(Seurat)

sref <- readRDS("../xCell2.0/reference_data/humanCellTypeAtlas.rds")

# all_ds <- unique(unlist(strsplit(unique(sref$series), ",")))
# all_ont <- unique(sref$cl.id)
all_label <- unique(sref$lowest.label)


# Endothelial cells - Epithelial cells
samples2use1 <- grep("endothe", sref$lowest.label)
labels.low1 <- sref$lowest.label[samples2use1]
unique(labels.low1)
labels.high1 <- rep("Endothelial cells", length(labels.low1))

samples2use2 <- grep("epi", sref$lowest.label)
labels.low2 <- sref$lowest.label[samples2use2]
unique(labels.low2)
labels.high2 <- rep("Epithelial cells", length(labels.low2))


data_sub <- sref$data[,c(samples2use1, samples2use2)]
dim(data_sub)
rm(sref)

sref.seurat <- CreateSeuratObject(counts = data_sub, project = "Super Ref")
sref.seurat <- FindVariableFeatures(sref.seurat, selection.method = "vst", nfeatures = 7000)
VariableFeaturePlot(sref.seurat)
sref.seurat <- ScaleData(sref.seurat, features = rownames(data_sub))

sref.seurat <- RunPCA(sref.seurat)
ElbowPlot(sref.seurat, ndims=45)

sref.seurat <- RunUMAP(sref.seurat, dims = 1:35)
sref.seurat <- AddMetaData(sref.seurat, metadata=labels, col.name = "labels")
DimPlot(sref.seurat, reduction = "umap", group.by = "labels", label = T)


