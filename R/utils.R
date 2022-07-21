library(tidyverse)
# This function return a pseudo-bulk expression matrix give a single cell RNA-Seq data

sc2pseudoBulk <- function(ref, labels, is_10x){
  if (class(ref) == "Seurat") {
    celltypes <- unique(labels$label)
    mat <- Seurat::GetAssayData(ref, assay = 'RNA', slot = 'data')
    for (ct in celltypes) {
      ct_mat <- mat[, labels$label == ct]
      if (!is_10x) {
        ct_mat[rowMeans(log2(ct_mat)) < 0.75,] <- 0 # Only for full length scRNA-Seq methods
      }
      if (ncol(ct_mat) >= 3) { # Minimum 3 cells per cluster to be a reference
        for (i in 1:5) { # Make five pseudo-bulk samples
          tmp_bulk <- rowSums(ct_mat[,sample(ncol(ct_mat), round(ncol(ct_mat)*0.5))]) # Use 50% of cells in the cluster each time
          if (i == 1L) {
            ct_mat_bulk <- tmp_bulk
          }else{
            ct_mat_bulk <- cbind(ct_mat_bulk, tmp_bulk)
          }
        }
        colnames(ct_mat_bulk) <- paste0(ct, ".", 1:5)
        if(ct==celltypes[1]){
          mat_bulk <- ct_mat_bulk
          lables_new <- data.frame("ont" = rep(unique(labels[labels$label == ct,1]), 5), "label" = rep(ct, 5), row.names = colnames(ct_mat_bulk))
        }else{
          mat_bulk <- cbind(mat_bulk, ct_mat_bulk)
          lables_new <- rbind(lables_new, data.frame("ont" = rep(unique(labels[labels$label == ct,1]), 5), "label" = rep(ct, 5), row.names = colnames(ct_mat_bulk)))
        }
      }else{
        print(paste0("WARNING: minimum 3 cells for ", ct, " required"))
      }
    }
    mat_bulk <- mat_bulk[rowSums(mat_bulk) > 0,]
  }else{
    # TODO: Make pseudo bulk for non-Seurat object
  }
  return(list("pseudoBulk" = as.matrix(mat_bulk), "newLabels" = lables_new))
}


# This function return a correlation matrix given the counts and cell types
getCellTypeCorrelation <- function(ref, labels){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  # Calculate median expression for each cell type
  median_expression <- lapply(celltypes, function(x){
    if (sum(samples == x) == 1) {
      as.numeric(ref[,samples == x])
    }else{
      apply(ref[,samples == x], 1, function(x) median(x, na.rm = TRUE))
    }
  })
  names(median_expression) <- celltypes

  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(median_expression[celltype_i][[1]], median_expression[celltype_j][[1]])
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(median_expression[celltype_i][[1]], median_expression[celltype_j][[1]])
  }

  return(cor_mat)
}


# This function return a vector of cell type dependencies
getDependencies <- function(OBOfile, labels){

  celltypes <- unique(labels[,2])
  cl <- suppressWarnings(ontologyIndex::get_ontology(OBOfile))

  dep_list <- vector(mode = "list", length = length(celltypes))
  names(dep_list) <- celltypes

  for (type in celltypes) {
    # Get cell type ontology
    ont <- as.character(unique(labels[labels[,2] == type, 1]))
    # Find descendants
    descendants <- ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = T)
    # Find ancestors
    ancestors <- ontologyIndex::get_ancestors(cl, terms = ont)
    ancestors <- ancestors[ancestors != ont]
    # Use only ontologies from the reference
    dep_cells <- c(descendants, ancestors)
    dep_cells <- dep_cells[dep_cells %in% labels[,1]]
    # Go back to cell type labels
    dep_cells <- unique(labels[labels[,1] %in% dep_cells, 2])

    dep_list[[type]] <- dep_cells
  }

  return(dep_list)
}


createInSilicoMixture <- function(ref, labels, cor_mat, mixture_fractions, add_noise){

  celltypes <- unique(labels[,2])

  fractions_mat_list <- lapply(celltypes, function(type) {

    # Get control vector
    control <- names(which.min(cor_mat[type,]))
    control_samples <- labels[,2] == control
    if (sum(control_samples) == 1) {
      control_vec <- as.vector(ref[,control_samples])
    }else{
      control_vec <- as.vector(apply(ref[,control_samples], 1, median))
    }

    # Get CTOI vector
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- as.vector(apply(ref[,type_samples], 1, median))
    }

    # Calculate fractions
    if (add_noise) {
      fractions_mat <- sapply(mixture_fractions, function(f) {
      (type_vec*f + control_vec*(1-f))*runif(1, 0.95, 1.05) # Multiply by a random number between 0.95-1.05 for noise
      })
    }else{
      fractions_mat <- sapply(mixture_fractions, function(f) {
      type_vec*f + control_vec*(1-f)
      })
    }
    colnames(fractions_mat) <- paste(mixture_fractions, control, sep = "_")
    fractions_mat
  })
  names(fractions_mat_list) <- celltypes

  # Transform list of matrices to one matrix
  ref_insilico_mat <- as.matrix(do.call(cbind.data.frame, fractions_mat_list))
  colnames(ref_insilico_mat) <- sub("\\.", "_", colnames(ref_insilico_mat))
  rownames(ref_insilico_mat) <- rownames(ref)

  # Split matrices into pure and fractions cell types
  pure_mat <- ref_insilico_mat[,grepl("_1_", colnames(ref_insilico_mat))]
  frac_mat <- ref_insilico_mat[,!grepl("_1_", colnames(ref_insilico_mat))]

  if (length(celltypes) != ncol(pure_mat)) {
    print("WARNING: Please remove _1_ from cell types names")
  }

  insilico_mat <- list("pureMat" = pure_mat, "fracMat" = frac_mat)

  return(insilico_mat)
}

# Generate a list with quantiles matrices for each cell type
makeQuantiles <- function(ref, labels, probs){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  quantiles_matrix <- lapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    # If there is one sample for this cell type - duplicate the sample to make a data frame
    if (sum(type_samples) == 1) {
      type.df <- cbind(ref[,type_samples], ref[,type_samples])
    }else{
      type.df <- ref[,type_samples]
    }
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, unique(c(probs, rev(1-probs))), na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes

  return(quantiles_matrix)
}


# Make score_mat tidy
# TODO: This function take too long
makeScoreMatTidy <- function(scores_mat){
  scores_mat_tidy <- scores_mat %>%
  as_tibble(., rownames = NA) %>%
  rownames_to_column(var = "signature") %>%
  pivot_longer(cols = -signature,
               names_to = c("mixture_ct", "mixture_fraction", "mixture_control"),
               names_sep = "_",
               values_to = "score") %>%
  separate(signature, into = "signature_ct", sep = "_", remove = FALSE, extra = "drop") # Speed bottleneck!!
}


# Plot signatures heatmap for a cell type
plotHeatMap <- function(type, scores_mat_tidy, signatures_collection_filtered = NULL, cor_mat){

  sig_score.df <- scores_mat_tidy %>%
    filter(signature_ct == type) %>%
    dplyr::select(signature, mixture_ct, score) %>%
    pivot_wider(names_from = mixture_ct, values_from = score) %>%
    column_to_rownames(var = "signature")

  sigs_to_use <- rownames(sig_score.df)
  celltype_order <- names(sort(cor_mat[type,], decreasing = TRUE))
  if (!is.null(signatures_collection_filtered)) {
    sigs_to_use <- sigs_to_use[sigs_to_use %in% names(signatures_collection_filtered)]
  }
  sig_score.df <- sig_score.df[sigs_to_use, celltype_order]
  sig_score.df <- sig_score.df[ ,colSums(is.na(sig_score.df)) == 0]

  pheatmap::pheatmap(sig_score.df, cluster_rows=F, cluster_cols=F, scale = "row", col= RColorBrewer::brewer.pal(11, "RdBu"))

}


# Helper function (for development) ------

# This function convents signatures_list from a nested lists to a list of GeneSet objects
makeGeneSetObjects <- function(signatures_list){
  all_signatures <- lapply(rapply(signatures_list, enquote, how="unlist"), eval)
  for (i in 1:length(all_signatures)) {
    all_signatures[[i]] <- GSEABase::GeneSet(all_signatures[i][[1]], setName = names(all_signatures[i]))
  }
  return(all_signatures)
}



# # Plot signatures heatmap for a cell type
# plotHeatMap <- function(ctoi, scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(1, 1, 1, 1), take_top_per = NULL, label = "ont"){
#
#   chnageLabel <- function(ref, ont, label = "fine"){
#     if (label == "main") {
#       out <- unique(ref$label.main[ref$label.ont == ont])
#       return(out[!is.na(out)][1])
#     }else if(label == "fine"){
#       out <- unique(ref$label.fine[ref$label.ont == ont])
#       return(out[!is.na(out)][1])
#     }
#   }
#
#   # Take only scores with 100% CTOI mixture
#   scores_mat_sub <- scores_mat[startsWith(rownames(scores_mat), ctoi), grep("*_1_", colnames(scores_mat))]
#   colnames(scores_mat_sub) <- unlist(lapply(strsplit(colnames(scores_mat_sub), "_"), "[", 1))
#
#   # Remove dep celltypes
#   types2use <- names(sort(celltype_cor_mat[ctoi,], decreasing = T))
#   types2use <- types2use[!types2use %in% dep_list[[ctoi]]]
#
#   # Sort signatures by rank
#   sigs_sorted <- signatures_ranked %>%
#     filter(signature_ct == ctoi) %>%
#     mutate(CT_rank = coalesce(CT_rank, CTOI_rank), diff_rank = coalesce(diff_rank, CTOI_rank)) %>%
#     mutate(final_rank = ranks_weights[1]*CTOI_rank + ranks_weights[2]*CT_rank + ranks_weights[3]*diff_rank + ranks_weights[4]*grubbs_rank) %>%
#     group_by(signature_ct) %>%
#     arrange(-final_rank, .by_group = TRUE)
#
#   # Filter signatures
#   if (!is.null(take_top_per)) {
#     sigs_sorted <- sigs_sorted %>%
#       filter(final_rank > quantile(final_rank, 1-take_top_per))
#   }
#
#   scores_mat_sub <- scores_mat_sub[sigs_sorted$signature, types2use]
#
#   # Change celltype labels
#   if (label == "main") {
#     colnames(scores_mat_sub) <- unlist(lapply(colnames(scores_mat_sub), function(x) chnageLabel(ref, x, label = "main")))
#   }else if(label == "fine"){
#     colnames(scores_mat_sub) <- unlist(lapply(colnames(scores_mat_sub), function(x) chnageLabel(ref, x)))
#   }
#
#   hm <- pheatmap::pheatmap(scores_mat_sub, cluster_rows=F, cluster_cols=F, scale = "row", col=brewer.pal(11, "RdBu")) # Change the yellow color
#   return(hm)
# }


# Ge DE genes with limma
getDEG <- function(ref, ctoi, ct2, plot = TRUE){
  library(limma)

  coldata.df <- ref@colData[ref@colData$label.ont %in% c(ctoi, ct2),]
  coldata.df$is_ctoi <- as.factor(ifelse(coldata.df$label.ont == ctoi, "yes", "no"))

  counts_mat <- ref@assays@data$logcounts[,rownames(coldata.df)]
  counts_mat <- counts_mat[rowSums(counts_mat) != 0,]

  design <- model.matrix(~is_ctoi, coldata.df)
  fit <- lmFit(counts_mat, design)
  fit <- eBayes(fit)

  toptable <- topTable(fit, n = Inf)
  toptable$Gene.symbol <- rownames(toptable)

  if (plot) {
    EnhancedVolcano::EnhancedVolcano(toptable,
                                     lab = toptable$Gene.symbol,
                                     x = 'logFC',
                                     y = 'P.Value')
  }


return(toptable)
}


# Change label

changeLabel <- function(ref, labels, to = "fine"){
  labels.df <- drop_na(as_tibble(colData(ref)))

  for (i in 1:ncol(labels.df)) {
    if(all(labels %in% labels.df[,i])){
      break
    }
  }

  if (to == "main") {
    out <- unname(unlist(lapply(labels, function(x) {unique(labels.df[labels.df[,i] == x, 1])})))
  }else if(to == "fine"){
    out <- unname(unlist(lapply(labels, function(x) {unique(labels.df[labels.df[,i] == x, 2])})))
    }else if(to == "ont"){
      out <- unname(unlist(lapply(labels, function(x) {unique(labels.df[labels.df[,i] == x, 3])})))
    }

  return(out)

}






