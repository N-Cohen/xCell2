# Main functions for xCell2

# TODO:

########################################################################################
# Generate cell-type lineage automatically
########################################################################################

# USAGE:
# ontology_table - a table with cell-type <label> column and ontology ID <ont>
#   For example:
# A tibble: 25 × 2
#      label                ont
#      <chr>               <chr>
#   1 CD4_T_cells         CL:0000624
#   2 CD8_T_cells         CL:0000625
#   3 T_helpers           CL:0000912
#   4 NK_cells            CL:0000623
# out_file - path to cell-type lineage file (.tsv)

# DEPENDENCIES:
# tidyverse, ontoProc, ontologyIndex

xCell2GetLineage <- function(ontology_table, out_file){

  cl <- ontoProc::getCellOnto()

  ontology_table %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(descendants = list(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)),
           ancestors = list(ontologyIndex::get_ancestors(cl, terms = ont))) %>%
    mutate(ancestors = list(ancestors[ancestors != ont])) %>%
    mutate(descendants = paste(pull(ct_ontologies[pull(ct_ontologies[,1]) %in% descendants, 2]), collapse = ", "),
           ancestors = paste(pull(ct_ontologies[pull(ct_ontologies[,1]) %in% ancestors, 2]), collapse = ", ")) %>%
    write_tsv(out_file)

  warning("It is recommended that you manually check the cell-type ontology file: ", out_file)
}



########################################################################################
# Train signatures
########################################################################################

# NOTE (!!!):
#  - No underscore in cell types labels
#  - Cell type onthology with colon in example: CL:0000545 (not CL_0000545)

# UPDATES:
# - in createSignatures I made cor_cells_cutoff to be the top 10% value of all correlation values
# - per_cells have been removed and is now 0.5 !!!
# - cor_mat is now in spearman (not pearson)

# DEPENDENCIES:
# tidyverse, pheatmap, singscore, GSEABase, outliers, ontoProc, ontologyIndex

# USAGE:
# ref -  ref matrix of the new reference (genes x samples) (!!!)
# (!!!) Check that the numbers are not strings.
# labels - a data frame with rows correspond to samples in ref:
#   (1) first column for cell type onthology
#   (2) second column for cell type name (charaters) (!!!)
#   (3) third column if samples should be in test (boolean)
# (!!!) If a cell-type doesn't have an ontology -> write the ontology of the most recent ancestor
# OBOfile - https://obofoundry.org/ontology/cl.html

# Remove
if (1 == 0) {
  data_type = "rnaseq";  score_method = "singscore"; mixture_fractions = c(0.001, seq(0.01, 0.25, 0.02), 1)
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  min_genes = 5; max_genes = 200; is_10x = TRUE
}


xCell2Train <- function(ref, labels, ontology_file_checked, data_type, score_method = "singscore", mixture_fractions = c(.001, seq(.01, .25, .02)),
                         probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                         min_genes = 7, max_genes = 200, is_10x = TRUE){


  # Validate inputs
  if (!"matrix" %in% class(ref)) {
    stop("ref should be as matrix!")
  }

  if (!"data.frame" %in% class(labels)) {
    stop("labels should be as dataframe!")
  }

  if (!data_type %in% c("rnaseq", "array", "scrnaseq")) {
    stop("data_type should be rnaseq, array or scrnaseq!")
  }

  underscore_check <- grepl("_", labels$label)    # If cell-types labels have underscore -> replace to dash
  if (sum(underscore_check) != 0) {
    message("Changing underscores to dashes in cell-types labels!")
    labels$label[underscore_check] <- gsub("_", "-", labels$label[underscore_check])
  }


  source("R/utils.R")
  # train/test samples
  test <- labels$is_test
  train <- !labels$is_test

  # Make pseudo-bulk for single-cell data
  if (data_type == "sc") {
    message("Making pseudo-bulk reference from scRNA-Seq data...")
    out <- sc2pseudoBulk(ref[,train], labels[train,], is_10x = is_10x)
    ref <- out$pseudoBulk
    labels <- out$newLabels
  }

  # (1) Make a table with median expression of pure cell types
  message("Calculating cell types median expression...")
  pure_ct_mat_train <- makePureCTMat(ref[,train], labels[train,])
  pure_ct_mat_test <- makePureCTMat(ref[,test], labels[test,])

  # (2) Build cell types correlation matrix
  message("Calculating cell type correlation matrix...")
  cor_mat <- getCellTypeCorrelation(pure_ct_mat_train)

  # (3) Get cell type dependencies list
  message("Finding cell types dependencies...")
  dep_list <- getDependencies(ontology_file_checked)

  # (4) Generate a list of quantiles matrices
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref[,train], labels[train,], probs)

  source("R/create_signatures.R")
  # (5) Generate signatures for each cell type
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref[,train], labels[train,], dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)

  source("R/filter_signatures.R")
  # (6) Filter signatures
  message("Filtering signatures...")
  filter_signature_out <- filterSignatures(pure_ct_mat_test, dep_list, signatures_collection, score_method, take_top_per = 0.1, max_sigs = 10)
  scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  signatures_collection_filtered <- filter_signature_out$sigCollectionFilt
  # plotHeatMap("CD8-T-cells", scores_mat_pure_tidy, signatures_collection_filtered = NULL, cor_mat)
  # plotHeatMap("CD8-T-cells", scores_mat_pure_tidy, signatures_collection_filtered = signatures_collection_filtered, cor_mat)

  # (7) Weight signatures with Elastic Net
  # TODO: change script name
  # source("train_models.R")
  source("R/train_models_en.R")
  models <- trainModels(ref[,test], labels[test,], pure_ct_mat_test, signatures_collection_filtered, dep_list, mixture_fractions, n_samples = 100, n_random_samples = 10)

  # source("transform_scores.R")
  # # (6) Fit signatures with Lasso regression
  # transformScores_out <- transformScores(ref, labels, signatures_collection_filtered, dep_list, cor_mat, mixture_fractions, score_method, add_noise = FALSE)
  # calibration_values <- transformScores_out[[1]]
  # scores_transformed <- transformScores_out[[2]]

  # Create S4 object for the new reference
  setClass("xCell2 Reference", slots = list(
    ref = "matrix",
    labels = "data.frame",
    correlationMatrix = "matrix",
    dependencies = "list",
    #scoresMatrix = "data.frame",
    signatures = "GeneSetCollection",
    models = "tbl"
    # calibrationValues = "data.frame"
  ))


  xCell2Ref.s4 <- new("xCell2 Reference", ref = ref, labels = labels, correlationMatrix = cor_mat, dependencies = dep_list, #scoresMatrix = as.data.frame(scores_mat_pure_tidy),
                      signatures = signatures_collection_filtered, models = models[,c(1,2,4)])

  return(xCell2Ref.s4)

}


# xCell2Analysis ---------------



xCell2Analysis <- function(mix, ref){

  celltypes <- unique(ref@labels$label)

  # Rank mixture genes
  mix_ranked <- singscore::rankGenes(mix)

  # # Score mixture with signatures
  # scores_list <- lapply(celltypes, function(ct){
  #   ct_sigs <- ref@signatures[startsWith(names(ref@signatures), ct)]
  #
  #   scores_out <- sapply(ct_sigs, function(sig){
  #     singscore::simpleScore(mix_ranked, upSet=sig, centerScore = FALSE)$TotalScore
  #   })
  # })
  # names(scores_list) <- celltypes

  # scores_list <- lapply(celltypes, function(ct){
  #   ct_sig_set <- pull(ref@model[ref@model$celltype == ct, 2])[[1]]
  #
  #   set_score_mat <- sapply(ct_sig_set, simplify = TRUE, function(set){
  #     singscore::simpleScore(mix_ranked, upSet = set, centerScore = FALSE)$TotalScore
  #   })
  # })


  # scores_list <- lapply(celltypes, function(ct){
  #   ct_sigs <- pull(ref@models[ref@models$celltype == ct,2])[[1]]
  #
  #       scores_out <- sapply(ct_sigs, function(sig){
  #     singscore::simpleScore(mix_ranked, upSet=sig, centerScore = FALSE)$TotalScore
  #   })
  # })
  # names(scores_list) <- celltypes


  # predictions <- ref@models %>%
  #   rowwise() %>%
  #   mutate(predictions = list(rowMeans(sapply(signatures, simplify = TRUE, function(set){
  #     singscore::simpleScore(mix_ranked, upSet = set, centerScore = FALSE)$TotalScore
  #   }))))
  #
  # xCell2_out <- matrix(unlist(predictions$predictions), ncol=ncol(mix), byrow=TRUE, dimnames = list(celltypes, colnames(mix)))



  scores_out <- sapply(ref@signatures, function(x){
    singscore::simpleScore(mix_ranked, upSet=x, centerScore = FALSE)$TotalScore
  })
  colnames(scores_out) <- names(ref@signatures)
  rownames(scores_out) <- colnames(mix)


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
