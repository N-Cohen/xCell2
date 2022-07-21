# Main functions for xCell2

# xCell2NewData ---------------

# TODO:
# Add function for xCell2NewData for single-cell RNA-seq reference
# Add slot for scores in the S4 object and also a slot for heatmaps with a boolean argument

# NOTE (!!!):
#  - No underscore in cell types labels
#  - Cell type onthology with colon in example: CL:0000545 (not CL_0000545)

# UPDATES:
# - in createSignatures I made cor_cells_cutoff to be the top 10% value of all correlation values
# - per_cells have been removed and is now 0.5 !!!

# USAGE:
# ref -  ref matrix of the new reference (genes x samples)
# labels - a data frame with row names as samples in ref, first column for cell type onthology and second column for cell type name (charaters) (!!!)
# OBOfile - https://obofoundry.org/ontology/cl.html


xCell2NewData <- function(ref, labels, data_type, mixture_fractions = c(.001, seq(.01, .25, .02), 1),
                         cor_cells_cutoff = .92, probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                         min_genes = 7, max_genes = 200, is_10x = TRUE){


  OBOfile <- "/home/almogangel/xCell2/data/cl.obo"
  setwd("/home/almogangel/xCell2/R/") # Remove

  source("utils.R")

  # Make pseudo-bulk for single-cell data
  if (data_type == "sc") {
    message("Making pseudo-bulk reference from scRNA-Seq data...")
    out <- sc2pseudoBulk(ref, labels, is_10x = is_10x)
    ref <- out$pseudoBulk
    labels <- out$newLabels
  }

  # If ref is a df -> convert to matrix
  if (!"matrix" %in% class(ref)) {
    ref <- as.matrix(ref)
  }

  # (1) Build cell types correlation matrix
  message("Calculating cell type correlation matrix...")
  cor_mat <- getCellTypeCorrelation(ref, labels)
  # saveRDS(cor_mat, "../data_for_dev/cor_mat.RData")
  # cor_mat <- readRDS("../data_for_dev/cor_mat.RData")

  # (2) Get cell type dependencies list
  message("Finding cell types dependencies...")
  dep_list <- getDependencies(OBOfile, labels)
  # saveRDS(dep_list, "../data_for_dev/dep_list.RData")
  # dep_list <- readRDS("../data_for_dev/dep_list.RData")

  # (3) Generate a list of quantiles matrices
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs)
  # saveRDS(quantiles_matrix, "../data_for_dev/quantiles_matrix.RData")
  # quantiles_matrix <- readRDS("../data_for_dev/quantiles_matrix.RData")

  source("create_signatures.R")
  # (4) Generate signatures for each cell type
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, cor_cells_cutoff, diff_vals, min_genes, max_genes)
  # saveRDS(signatures_collection, "../data_for_dev/signatures_collection.RData")
  # signatures_collection <- readRDS("../data_for_dev/signatures_collection.RData")

  # (5) Build in silico mixtures
  message("Building in silico mixtures...")
  inSilico_mat <- createInSilicoMixture(ref, labels, cor_mat, mixture_fractions, add_noise = FALSE)
  # saveRDS(inSilico_mat, "../data_for_dev/inSilico_mat.RData")
  # saveRDS(inSilico_mat, "../data_for_dev/inSilico_mat_withNoise.RData")
  # inSilico_mat <- readRDS("../data_for_dev/inSilico_mat.RData")

  source("score_signatures.R")
  source("filter_signatures.R")
  # (6) Filter signatures
  message("Filtering signatures...")
  scores_mat_pure <- scoreSignatures(signatures_collection, inSilico_mat$pureMat, dep_list, score_method = "singscore")
  scores_mat_pure_tidy <- makeScoreMatTidy(scores_mat_pure)
  signatures_collection_filtered <- filterSignatures(signatures_collection, scores_mat_pure_tidy, take_top_per = 0.1)
  # saveRDS(signatures_filtered, "signatures_filtered.RData")
  # signatures_filtered <- readRDS("signatures_filtered.RData")
  # plotHeatMap(labels$label[1], scores_mat_pure_tidy, mixture_frac = 1, filter_signatures = NULL, cor_mat)
  # plotHeatMap(labels$label[1], scores_mat_pure_tidy, mixture_frac = 1, filter_signatures = signatures_collection_filtered, cor_mat)

  source("transform_scores.R")
  # (7) Fit signatures with Lasso regression
  scores_mat_frac <- scoreSignatures(signatures_collection_filtered, inSilico_mat$fracMat, dep_list, score_method = "singscore")
  scores_mat_frac_tidy <- makeScoreMatTidy(scores_mat_frac)
  transformScores_out <- transformScores(scores_mat_frac_tidy)
  calibration_values <- transformScores_out[[1]]
  scores_transformed <- transformScores_out[[2]]
  # saveRDS(calibration_values, "calibration_values.RData")
  # saveRDS(scores_transformed, "scores_transformed.RData")
  # calibration_values <- readRDS("calibration_values.RData")
  # scores_transformed <- readRDS("scores_transformed.RData")


  # Create S4 object for the new reference
  setClass("xCell2 Reference", slots = list(
    ref = "matrix",
    labels = "data.frame",
    correlationMatrix = "matrix",
    dependencies = "list",
    scoresMatrix = "data.frame",
    signatures = "GeneSetCollection",
    calibrationValues = "data.frame"
  ))


  xCell2Ref.s4 <- new("xCell2 Reference", ref = ref, labels = labels, correlationMatrix = cor_mat, dependencies = dep_list, scoresMatrix = as.data.frame(scores_mat_pure_tidy),
                      signatures = signatures_collection_filtered, calibrationValues = data.frame(calibration_values))
  # saveRDS(xCell2Ref.s4, "../data_for_dev/BlueprintEncode/xCell2RefBlueprintEncode.RData")

  return(xCell2Ref.s4)

}


# xCell2Analysis ---------------



xCell2Analysis <- function(mix, ref){

  # Rank mixture genes
  mix_ranked <- singscore::rankGenes(mix)

  # Score mixture with signatures
  scores_out <- sapply(ref@signatures, function(x){
    singscore::simpleScore(mix_ranked, upSet=x, centerScore = FALSE)$TotalScore
  })
  colnames(scores_out) <- names(ref@signatures)
  rownames(scores_out) <- colnames(mix)

  # Make results tidy
  xCell2_out <- scores_out %>%
    as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to	= c("signature_ct", "quantile", "diff", "n_sim", "n_passed"), values_to = "score", names_sep = "_") %>%
    group_by(Sample, signature_ct) %>%
    summarise(mean_score = mean(score)) %>%
    left_join(ref@calibrationValues, by = "signature_ct") %>%
    mutate(adjusted_score = (mean_score-shift_value)*scale_factor) %>%
    dplyr::select(Sample, signature_ct, adjusted_score) %>%
    pivot_wider(names_from = Sample, values_from = adjusted_score) %>%
    data.frame(row.names = 1, check.names = FALSE) %>%
    as.matrix()

  return(xCell2_out)
}
