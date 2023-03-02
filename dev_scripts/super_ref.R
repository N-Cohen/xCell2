library(tidyverse)

mref <-  readRDS("/bigdata/mahmoudy/humanCellTypeAtlas.rds")
aaaa
dim(mref$data)
colnames(mref$data)[1:10]
View(mref$data[1:10, 1:10])
length(mref$tissue_titles)
unique(mref$lowest.label)

all_ds <- unique(unlist(strsplit(unique(mref$series), ",")))
unique(all_models_annot$Dataset)[!unique(all_models_annot$Dataset) %in% all_ds]


ds2check <- all_ds[! all_ds %in% unique(all_models_annot$Dataset)]
mref$cl.id[mref$series == ds2check[1]]


cl <- ontoProc::getCellOnto()
related.cl <- c(ontologyIndex::get_descendants(cl, roots = "CL:0002543", exclude_roots = TRUE),
  ontologyIndex::get_ancestors(cl, terms = "CL:0002543"))
related.cl <- unique(mref$cl.id[mref$cl.id %in% related.cl])

mref$tissue_titles[mref$cl.id %in% related.cl]

data_sub <- mref$data[,mref$cl.id %in% related.cl]
library(Seurat)
pca <- RunPCA(nsclc.mLN.small)
