library(tidyverse)
library(RColorBrewer)
library(dichromat)

setwd("data_for_dev/BlueprintEncode/")


# Input:
ref = readRDS("BE.RData")
OBOfile = "~/Documents/xCell2.0/cl.obo"
use_genes = NULL
RNA_seq = TRUE
labels = "ont"
mixture_fractions = c(seq(0.001, 0.009, 0.001), seq(.01,.25, .01), 1)


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
# BlueprintEncode already in log counts !!! - https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf
# if (RNA_seq) {
#   counts <- log2(counts + 1)
#   message("Data type: RNA-Seq (using log2 transformation for counts)")
# }else{
#   message("Data type: array")
# }


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


source("../../R/utils.R")
# (5) Get correlation matrix between cell types
# celltype_cor_mat <- getCellTypeCorrelation(counts, samples, celltypes)
# saveRDS(celltype_cor_mat, "celltype_cor_mat.RData")
celltype_cor_mat <- readRDS("celltype_cor_mat.RData")


# (6) Get cell type dependencies list
# dep_list <- getDependencies(OBOfile, celltypes)
# saveRDS(dep_list, "dep_list.RData")
dep_list <- readRDS("dep_list.RData")


# (7) Get in silico mixtures
# inSilico_mat <- createInSilicoMixture(counts, celltypes, celltype_cor_mat, fractions = mixture_fractions)
# saveRDS(inSilico_mat, "inSilico_mat.RData")
inSilico_mat <- readRDS("inSilico_mat.RData")


source("../../R/create_signatures.R")
# (8) Get a list of signatures for each cell type
# signatures_list <- createSignatures(counts, samples, celltypes, dep_list = dep_list, celltype_cor_mat = celltype_cor_mat, cor_cells_cutoff = .92)
# saveRDS(signatures_list, "signatures_list.RData")
signatures_list <- readRDS("signatures_list.RData")


source("../../R/score_signatures.R")
# (9) Score all signatures vs. the in silico matrix
# all_signatures <- makeGeneSetObjects(signatures_list)
# saveRDS(all_signatures, "all_signatures.RData")
all_signatures <- readRDS("all_signatures.RData")
# scores_mat <- scoreSignatures(all_signatures, celltypes, inSilico_mat, dep_list)
# saveRDS(scores_mat, "scores_mat.RData")
scores_mat <- readRDS("scores_mat.RData")


source("../../R/rank_signatures.R")
# (10) Rank signatures
# signatures_ranked <- rankSignatures(scores_mat)
# saveRDS(signatures_ranked, "signatures_ranked.RData")
signatures_ranked <- readRDS("signatures_ranked.RData")


# # Example heatmaps (Remove)
# plotHeatMap("CL:0000136", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(0, 0, 0, 1), label = "fine")
# plotHeatMap("CL:0000863", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(0, 0, 0, 1), label = "fine")
# plotHeatMap("CL:0000863", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(0, 0, 0, 0), label = "fine")
# plotHeatMap("CL:0000624", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(0, 0, 0, 0), label = "fine")
# plotHeatMap("CL:0000624", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(.5, .4, .1), take_top_per = .1, label = "fine")
# plotHeatMap("CL:0000786", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(.2, .4, .5), label = "fine")
# plotHeatMap("CL:0000786", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(.5, .4, .1), take_top_per = .1, label = "fine")
# # Something wrong with Fibroblasts
# signatures_ranked %>%
#   filter(signature_ct == "CL:0000057") %>%
#   arrange(-grubbs_rank)
# plotHeatMap("CL:0000057", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(0, 0, 0, 0), label = "fine")
# plotHeatMap("CL:0000057", scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(.5, .4, .1), take_top_per = .1, label = "fine")


source("../../R/filter_signatures.R")
# (11) Filter signatures

# TODO: If you use only grubbs_rank to filter signatures then you might lose some signature because all the pvalues of those signatures are the same (0) (!!!)

# signatures_filtered <- filterSignatures(signatures_ranked, ranks_weights = c(.1, .1, .1, 1), take_top_per = .1)
# saveRDS(signatures_filtered, "signatures_filtered.RData")
signatures_filtered <- readRDS("signatures_filtered.RData")


source("../../R/fit_signatures.R")
# (12) Fit signatures with Lasso regression
# signatures_fit <- fitSignatures(signatures_filtered, scores_mat)
# saveRDS(signatures_fit, "signatures_fit.RData")
signatures_fit <- readRDS("signatures_fit.RData")



source("../../R/spillover_correction.R")
# (13) Build a spillover correction matrix
spill_mat <- spillOverMat()
# saveRDS(spill_mat, "spill_matRData")
# signatures_fit <- readRDS("spill_mat.RData")




xCell2trainRef <- function(ref, use_genes, OBOfile, RNA_seq, labels, mixtures_fractions, cor_cells_cutoff = .92, ranks_weights = c(.1, .1, .1, 1), take_top_per = .1){


  # Check if reference is in SummarizedExperiment object
  if(class(ref) != "SummarizedExperiment"){
    stop(paste0("Reference dataset is not in a SummarizedExperiment class."))
  }

  # Use only genes in "use_genes"
  if (!is.null(use_genes)) {
    counts <- ref@assays@data$logcount[use_genes,]
    message(paste0("Using a subset of ", length(use_genes), " genes."))
  }else{
    counts <- ref@assays@data$logcount
  }

  # Choose cell type annotation
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


  source("../../R/utils.R")
  message("Calculating correlation matrix...")
  celltype_cor_mat <- getCellTypeCorrelation(counts, samples, celltypes)
  dep_list <- getDependencies(OBOfile, celltypes)
  inSilico_mat <- createInSilicoMixture(counts, celltypes, celltype_cor_mat, fractions = mixtures_fractions)


  source("../../R/create_signatures.R")
  signatures_list <- createSignatures(counts, samples, celltypes, dep_list, celltype_cor_mat, cor_cells_cutoff = cor_cells_cutoff)


  source("../../R/score_signatures.R")
  scores_mat <- scoreSignatures(signatures_list, inSilico_mat, dep_list)


  source("../../R/rank_signatures.R")
  signatures_ranked <- rankSignatures(scores_mat)


  source("../../R/filter_signatures.R")
  signatures_filtered <- filterSignatures(signatures_ranked, ranks_weights = ranks_weights, take_top_per = take_top_per)


  source("../../R/fit_signatures.R")
  signatures_fit <- fitSignatures(signatures_filtered, scores_mat)


  source("../../R/spillover_correction.R")



  # Build S4 object
  setClass("xCell2TrainedReference",
           slots = list(counts = "matrix", samples = "character", celltypes = "character"))

  xcell2_ref <- new("xCell2TrainedReference", counts = counts, samples = samples, celltypes = celltypes)


  return(xcell2_ref)

}

