# Main functions for xCell2

# TODO:

########################################################################################
# Generate cell-type lineage automatically
########################################################################################

# USAGE:
# labels - a data frame with rows correspond to samples in ref:
#   (1) first column for cell type onthology
#   (2) second column for cell type name (charaters)

# DEPENDENCIES:
# tidyverse, ontoProc, ontologyIndex

xCell2GetLineage <- function(labels, out_file){

  cl <- ontoProc::getCellOnto()

  labels_uniq <- labels %>%
    as_tibble() %>%
    unique()

  labels_uniq %>%
    rowwise() %>%
    mutate(descendants = list(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)),
           ancestors = list(ontologyIndex::get_ancestors(cl, terms = ont))) %>%
    mutate(ancestors = list(ancestors[ancestors != ont])) %>%
    mutate(descendants = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% descendants, 2]), collapse = ", "),
           ancestors = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% ancestors, 2]), collapse = ", ")) %>%
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
#   (2) second column for cell type name (charaters)


# Remove
if (1 == 0) {
  data_type = "rnaseq";  score_method = "singscore"; mixture_fractions = c(0.001, seq(0.01, 0.25, 0.02), 1)
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  min_genes = 5; max_genes = 500; is_10x = TRUE
}


xCell2Train <- function(ref, labels, ontology_file_checked, data_type, score_method = "singscore", mixture_fractions = c(.001, seq(.01, .25, .02)),
                         probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                         min_genes = 7, max_genes = 500, is_10x = TRUE){


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

  # Make pseudo-bulk for single-cell data
  if (data_type == "sc") {
    message("Making pseudo-bulk reference from scRNA-Seq data...")
    out <- sc2pseudoBulk(ref, labels, is_10x = is_10x)
    ref <- out$pseudoBulk
    labels <- out$newLabels
  }

  # (1) Make a table with median expression of pure cell types
  message("Calculating cell types median expression...")
  pure_ct_mat <- makePureCTMat(ref, labels)

  # (2) Build cell types correlation matrix
  message("Calculating cell type correlation matrix...")
  cor_mat <- getCellTypeCorrelation(pure_ct_mat)

  # (3) Get cell type dependencies list
  message("Finding cell types dependencies...")
  dep_list <- getDependencies(ontology_file_checked)

  # (4) Generate a list of quantiles matrices
  message("Calculating quantiles...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs)

  source("R/create_signatures.R")
  # (5) Generate signatures for each cell type
  message("Generating signatures...")
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)

  source("R/filter_signatures.R")
  # (6) Filter signatures
  message("Filtering signatures...")
  filter_signature_out <- filterSignatures(pure_ct_mat = pure_ct_mat, dep_list, signatures_collection, score_method, take_top_per = 0.1, max_sigs = 10)
  scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  signatures_collection_filtered <- filter_signature_out$sigCollectionFilt
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = NULL, cor_mat)
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = signatures_collection_filtered, cor_mat)

  # (7) Weight signatures with Elastic Net
  # TODO: change script name
  # source("train_models.R")
  #source("R/train_models_tmp.R")
  #models <- trainModels(ref, labels, dep_list, pure_ct_mat_test, signatures_collection_filtered, mixture_fractions)


  # Create S4 object for the new reference
  setClass("xCell2 Reference", slots = list(
    ref = "matrix",
    labels = "data.frame",
    correlationMatrix = "matrix",
    dependencies = "list",
    signatures = "GeneSetCollection"
    # models = "tbl"
  ))


  xCell2Ref.s4 <- new("xCell2 Reference", ref = ref, labels = labels, correlationMatrix = cor_mat, dependencies = dep_list,
                      signatures = signatures_collection_filtered)

  return(xCell2Ref.s4)

}



########################################################################################
# Analyze new mixtures
########################################################################################


xCell2Analysis <- function(mix, ref){

  scoreMixtures <- function(ctoi, mixture, signatures_ctoi){

    mixture_ranked <-  singscore::rankGenes(mixture)
    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })


    return(scores)
  }


  celltypes <- unique(ref@labels$label)

  # Rank mixture genes
  mix_ranked <- singscore::rankGenes(mix)

  xCell2_out <- ref@models %>%
    rowwise() %>%
    mutate(scores = list(scoreMixtures(ctoi = label, mixture = mix_ranked, signatures_ctoi = signatures))) %>%
    # For lasso:
    #mutate(predictions = list(as.numeric(predict(lasso, newx = scores)))) %>%
    # For GGRF:
    #mutate(predictions = list(as.numeric(predict(GRRF, newdata = scores)))) %>%
    dplyr::select(label, predictions) %>%
    unnest(predictions) %>%
    mutate(samples = rep(colnames(mix), length(celltypes))) %>%
    pivot_wider(names_from = samples, values_from = predictions) %>%
    as.data.frame()

  rownames(xCell2_out) <- xCell2_out[,1]
  xCell2_out <- xCell2_out[,-1]



  return(xCell2_out)
}
