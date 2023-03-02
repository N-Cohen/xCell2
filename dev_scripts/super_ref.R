library(tidyverse)
library(Seurat)

cl <- ontoProc::getCellOnto()

sref <- readRDS("/bigdata/mahmoudy/humanCellTypeAtlas.rds")


dim(sref$data)
colnames(sref$data)[1:10]
rownames(sref$data)[1:10]
View(sref$data[1:10, 1:10])
length(sref$tissue_titles)
unique(sref$lowest.label)

all_ds <- unique(unlist(strsplit(unique(sref$series), ",")))
all_ont <- unique(sref$cl.id)

ont <- all_ont[1]

ont2check <- c()
ont2check_old <- c()

while (length(ont2check) != length(ont2check_old)) {
  related.cl <- c(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE),
                  ontologyIndex::get_ancestors(cl, terms = ont))
  related.cl <- unique(sref$cl.id[sref$cl.id %in% related.cl])
  ont2check <- c(ont2check, related.cl)

}


unique(sref$lowest.label[sref$cl.id %in% related.cl])

data_sub <- sref$data[,sref$cl.id %in% related.cl]

sref.seurat <- CreateSeuratObject(counts = sref$data, project = "Super Ref")
sref.seurat <- FindVariableFeatures(sref.seurat, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(sref.seurat)
sref.seurat <- ScaleData(sref.seurat, features = rownames(data_sub))

# PCA
sref.seurat <- RunPCA(sref.seurat)
DimPlot(sref.seurat, reduction = "pca")
ElbowPlot(sref.seurat, ndims=35)

# UMAP
sref.seurat <- RunUMAP(sref.seurat, dims = 1:15)
sref.seurat <- AddMetaData(sref.seurat, metadata=sref$lowest.label[sref$cl.id %in% related.cl], col.name = "labels")
DimPlot(sref.seurat, reduction = "umap", group.by = "labels", label = T)


