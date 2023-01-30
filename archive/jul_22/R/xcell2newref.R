# xCell2NewRef function
library(GSEABase)
library(tidyverse)

setwd("~/Documents/xCell2.0/xCell2/R/")
ref <- readRDS("../data_for_dev/BlueprintEncode/old/BE.RData")

# Parameters:
# counts is a gene x samples counts matrix
counts = ref@assays@data$logcount
# labels is a df with first column ontologies, second column is the cell type label and row names are the samples of counts
labels = data.frame("ont" = ref$label.ont, "label" = ref$label.fine, row.names = colnames(ref))
# CL gene ontology OBO file
OBOfile = "../../cl.obo"
# Which genes to use in counts?
use_genes = NULL
# Data type of counts: "rnaseq", "array" or "scrnaseq"
data_type = "rnaseq"
# Fraction of cell type mixture
mixture_fractions = c(.001, seq(.01, .25, .01), 1)
# Add noise for each gene in the in silico mixtures
add_noise = FALSE
# Log 2 transform the counts? (only for RNA-Seq data)
log2transform = FALSE
# Consider cell types similar above this correlation value
cor_cells_cutoff = .92
# Probabilities for quantiles
probs = c(.1, .25, .333, .5)
# Difference between lower to upper quantiles
diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
# Minimum genes in a signature
min_genes = 7
# Maximum genes in a signature
max_genes = 200
# Minimum percent of cell types that pass the quantiles criteria for a gene to count as marker
per_cells = 0.95
# Define genes type of counts (for more details check GSEABase documentation)
genes_type = SymbolIdentifier()
# Define method to score signatures ("singscore"/"ssgsea")
score_method = "singscore"
# Filter top percentage signatures
take_top_per = 0.1


xCell2NewRef <- function(ref_name, counts, labels, OBOfile, use_genes, RNA_seq, mixture_fractions, log2transform,
                         cor_cells_cutoff = .92, probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                         min_genes = 7, max_genes = 200, per_cells = .95){

  # Use only genes in "use_genes"
  if (!is.null(use_genes)) {
    counts <- ref[use_genes,]
  }


  # log2 transformation for RNA-Seq data
  # Note: BlueprintEncode reference is already in log count - https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf
  if (RNA_seq & log2transform) {
    counts <- log2(counts + 1)
  }


  source("utils.R")
  # Build cell types correlation matrix
  message("Calculating cell type correlation matrix...")
  # cor_mat <- getCellTypeCorrelation(counts, labels)
  # saveRDS(cor_mat, "../data_for_dev/cor_mat.RData")
  cor_mat <- readRDS("../data_for_dev/cor_mat.RData")

  # Get cell type dependencies list
  message("Finding cell types dependencies...")
  # dep_list <- getDependencies(OBOfile, labels)
  # saveRDS(dep_list, "../data_for_dev/dep_list.RData")
  dep_list <- readRDS("../data_for_dev/dep_list.RData")

  # Build in silico mixtures
  message("Building in silico mixtures...")
  # inSilico_mat <- createInSilicoMixture(counts, labels, cor_mat, mixture_fractions, add_noise)
  # saveRDS(inSilico_mat, "../data_for_dev/inSilico_mat.RData")
  # saveRDS(inSilico_mat, "../data_for_dev/inSilico_mat_withNoise.RData")
  inSilico_mat <- readRDS("../data_for_dev/inSilico_mat.RData")

  # Generate a list of quantiles matrices
  message("Calculating quantiles...")
  # quantiles_matrix <- makeQuantiles(counts, labels, probs)
  # saveRDS(quantiles_matrix, "../data_for_dev/quantiles_matrix.RData")
  quantiles_matrix <- readRDS("../data_for_dev/quantiles_matrix.RData")



  source("create_signatures.R")
  # Generate signatures for each cell type
  message("Generating signatures...")
  # signatures_collection <- createSignatures(counts, labels, quantiles_matrix, probs, cor_mat, diff_vals, per_cells, genes_type, min_genes, max_genes)
  # saveRDS(signatures_collection, "../data_for_dev/signatures_collection.RData")
  signatures_collection <- readRDS("../data_for_dev/signatures_collection.RData")


  source("score_signatures.R")
  # Score all signatures vs. in silico mixtures
  message("Scoring mixtures...")
  # TODO: Fix warnings with ssGSEA
  # scores_mat <- scoreSignatures(signatures_collection, inSilico_mat, dep_list, score_method = "singscore")
  # saveRDS(scores_mat, "scores_mat_singscore.RData")
  # saveRDS(scores_mat, "scores_mat_ssgsea.RData")
  # scores_mat <- readRDS("scores_mat_singscore.RData")
  scores_mat <- readRDS("scores_mat_ssgsea.RData")

  # Make score_mat tidy
  # scores_mat_tidy <- makeScoreMatTidy(scores_mat)
  # saveRDS(scores_mat_tidy, "scores_mat_tidy.RData")
  # saveRDS(scores_mat_tidy, "scores_mat_tidy_ssgsea.RData")
  scores_mat_tidy <- readRDS("scores_mat_tidy.RData")


  # - Note: rank_signatures.R and filter_signatures.R are now merged in filter_signatures.R
  source("filter_signatures.R")
  # Filter signatures
  message("Filtering signatures...")
  # signatures_filtered <- filterSignatures(scores_mat_tidy, take_top_per)
  # saveRDS(signatures_filtered, "signatures_filtered.RData")
  signatures_filtered <- readRDS("signatures_filtered.RData")
  # plotHeatMap(labels$label[1], scores_mat_tidy, mixture_frac = 1, filter_signatures = NULL, cor_mat)
  # plotHeatMap(labels$label[1], scores_mat_tidy, mixture_frac = 1, filter_signatures = signatures_filtered, cor_mat)


  source("transform_scores.R")
  # Fit signatures with Lasso regression
  # transformScores_out <- transformScores(scores_mat_tidy, signatures_filtered)
  # calibration_values <- transformScores_out[[1]]
  # scores_transformed <- transformScores_out[[2]]
  # saveRDS(calibration_values, "calibration_values.RData")
  # saveRDS(scores_transformed, "scores_transformed.RData")
  calibration_values <- readRDS("calibration_values.RData")
  scores_transformed <- readRDS("scores_transformed.RData")


  # Create S4 object for the new reference
  setClass("xCell2Ref", slots =list(
    counts = "matrix",
    labels = "data.frame",
    correlationMatrix = "matrix",
    dependencies = "list",
    signatures = "GeneSetCollection",
    calibrationValues ="data.frame"
  ))


  xCell2Ref.s4 <- new("xCell2Ref", counts = counts, labels = labels, correlationMatrix = cor_mat, dependencies = dep_list,
                      signatures = signatures_collection[pull(signatures_filtered, signature)], calibrationValues = data.frame(calibration_values)[,c(1,3,5)])
  # saveRDS(xCell2Ref.s4, "../data_for_dev/BlueprintEncode/xCell2RefBlueprintEncode.RData")

  return(xCell2Ref.s4)

}
