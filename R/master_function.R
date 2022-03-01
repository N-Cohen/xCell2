source("R/utils.R")
source("R/create_signatures.R")

# Input:
ref = readRDS("~/Documents/xCell2.0/HPCA.RData")
OBOfile = "~/Documents/xCell2.0/cl.obo"
use_genes = NULL
RNA_seq = FALSE
labels = "ont"

# Master function steps:

# (1) Check if reference is in SummarizedExperiment object
if(class(ref) != "SummarizedExperiment"){
  stop(paste0("Reference dataset is not in a SummarizedExperiment class."))
}

# (2) Use only genes in "use_genes"
if (!is.null(use_genes)) {
  counts <- ref@assays@data$logcount[use_genes,]
}else{
  counts <- ref@assays@data$logcount
}


# (3) If reference is RNA-seq make log2 transformation
if (RNA_seq) {
  counts <- log2(counts + 1)
  message("Data type: RNA-Seq (using log2 transformation for counts)")
}else{
  message("Data type: array")
}


# (4) Train signatures for "main", "fine" or "CL ontology" cell types
# TODO: for now we only used "ont" - check the rest
if (labels == "ont") {
  samples <- ref$label.ont
  message("Cell type annotation: CL ontology")
}else if(labels == "main"){
  samples <- ref$label.main
  message("Cell type annotation: main")
}else if(labels == "fine"){
  samples <- ref$label.fine
  message("Cell type annotation: fine")
}
celltypes <- unique(samples[!is.na(samples)])


# (5) Get correlation matrix between cell types
# celltype_cor_mat <- getCellTypeCorrelation(counts, samples, celltypes)
celltype_cor_mat <- readRDS("data_for_dev/celltype_cor_mat.RData")


# (6) Get cell type dependencies list
# dep_list <- getDependencies(OBOfile, celltypes)
dep_list <- readRDS("data_for_dev/dep_list.RData")


# (7) Get in silico mixtures
# inSilico_mat <- createInSilicoMixture(counts, celltypes, celltype_cor_mat)
inSilico_mat <- readRDS("data_for_dev/inSilico_mat.RData")


# (8) Get a list of signatures for each cell type
# signatures_list <- createSignatures(counts, samples, celltypes, dep_list = dep_list, celltype_cor_mat = celltype_cor_mat, cor_cells_cutoff = .92)
signatures_list <- readRDS("data_for_dev/signatures_list.RData")


# (9) Score all signatures vs. the in silico matrix






