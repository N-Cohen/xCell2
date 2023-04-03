setwd(dir = "/Users/noam/Desktop/Technion/bioinformatics project/xCell2")

xCell2Train <- function(ref, labels, ontology_file_checked, data_type, mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02)),
                        probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                        min_genes = 5, max_genes = 500, is_10x = FALSE){
  
  
  # Validate inputs
  if (!any(class(ref) %in% c("matrix", "dgCMatrix", "Matrix"))) {
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
  
  # Build cell types correlation matrix
  message("Calculating cell-type correlation matrix...")
  pure_ct_mat <- makePureCTMat(ref, labels)
  cor_mat <- getCellTypeCorrelation(pure_ct_mat)
  
  # Get cell type dependencies list
  message("Loading dependencies...")
  dep_list <- getDependencies(ontology_file_checked)
  
  source("R/create_signatures.R")
  # Generate signatures for each cell type
  message("Generating signatures...")
  quantiles_matrix <- makeQuantiles(ref, labels, probs)
  signatures_collection <- createSignatures(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes)
  
  source("R/filter_signatures.R")
  # Filter signatures
  message("Filtering signatures...")
  filter_signature_out <- filterSignatures(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff = 0.8, simulations_cutoff = 0.8)
  scores_mat_pure_tidy <- filter_signature_out$scoreMatTidy
  signatures_collection_filtered <- filter_signature_out$sigCollectionFilt
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = NULL, cor_mat)
  # plotHeatMap("Neutrophils", scores_mat_pure_tidy, signatures_collection_filtered = signatures_collection_filtered, cor_mat)
  
  # TODO: (6) Weight signatures with Elastic Net
  # source("train_models.R")
  # source("R/train_models_tmp.R")
  # models <- trainModels(ref, labels, dep_list, pure_ct_mat_test, signatures_collection_filtered, mixture_fractions)
  
  
  # Create S4 object for the new reference
  setClass("xCell2 Reference", slots = list(
    ref = "matrix",
    labels = "data.frame",
    correlationMatrix = "matrix",
    dependencies = "list",
    score_mat = "tbl",
    all_signatures = "GeneSetCollection",
    filtered_signatures = "GeneSetCollection"
    # models = "tbl"
  ))
  
  
  xCell2Ref.S4 <- new("xCell2 Reference", ref = ref, labels = labels, correlationMatrix = cor_mat, dependencies = dep_list, score_mat = scores_mat_pure_tidy,
                      all_signatures = signatures_collection, filtered_signatures = signatures_collection_filtered)
  
  return(xCell2Ref.S4)
  
}
