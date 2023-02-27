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

# Generate a matrix of median expression of pure cell types
makePureCTMat <- function(ref, labels){

  celltypes <- unique(labels$label)

  pure_ct_mat<- sapply(celltypes, function(type){
    type_samples <- labels[,2] == type
    if (sum(type_samples) == 1) {
      type_vec <- as.vector(ref[,type_samples])
    }else{
      type_vec <- Rfast::rowMedians(ref[,type_samples])
    }
  })
  rownames(pure_ct_mat) <- rownames(ref)

  return(pure_ct_mat)
}

# This function return a correlation matrix given the counts and cell types
getCellTypeCorrelation <- function(pure_ct_mat){

  celltypes <- colnames(pure_ct_mat)

  # Use top 3K highly variable genes
  mean_gene_expression <- Rfast::rowmeans(pure_ct_mat)
  high_gene_expression_cutoff <- quantile(mean_gene_expression, 0.5, na.rm=TRUE) # Cutoff for top 50% expression genes
  top_expressed_gene <- mean_gene_expression > high_gene_expression_cutoff
  genes_sd <- apply(pure_ct_mat[top_expressed_gene,], 1, sd)
  sd_cutoff <-  sort(genes_sd, decreasing = TRUE)[1001]# Get top 1K genes with high SD as highly variable genes
  pure_ct_mat <- pure_ct_mat[genes_sd > sd_cutoff,]

  # Make correlation matrix
  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  lower_tri_coord <- which(lower.tri(cor_mat), arr.ind = TRUE)

  for (i in 1:nrow(lower_tri_coord)) {
    celltype_i <- rownames(cor_mat)[lower_tri_coord[i, 1]]
    celltype_j <- colnames(cor_mat)[lower_tri_coord[i, 2]]
    cor_mat[lower_tri_coord[i, 1], lower_tri_coord[i, 2]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
    cor_mat[lower_tri_coord[i, 2], lower_tri_coord[i, 1]] <- cor(pure_ct_mat[,celltype_i], pure_ct_mat[,celltype_j], method = "spearman")
  }

  return(cor_mat)
}


# This function return a vector of cell type dependencies
getDependencies <- function(ontology_file_checked){
  ont <- read_tsv(ontology_file_checked, show_col_types = FALSE)

  celltypes <- pull(ont[,2])
  dep_list <- vector(mode = "list", length = length(celltypes))
  names(dep_list) <- celltypes

  for (i in 1:nrow(ont)) {
    dep_cells <- c(strsplit(pull(ont[i,3]), ";")[[1]], strsplit(pull(ont[i,4]), ";")[[1]])
    dep_list[[i]] <- dep_cells[!is.na(dep_cells)]
  }

  return(dep_list)
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


makeMixture <- function(ctoi, ref, labels, pure_ct_mat, dep_list, n_samples){

  ctoi_mat <- ref[,labels[,2] == ctoi]

  # Make pure CTOI fractions matrix with n_samples (columns)
  pure_ctoi_mat <- sapply(1:n_samples, function(i){

    # Sample 0.2 of the #CTOI samples available in ref
    random_type_samples <- sample(1:ncol(ctoi_mat),
                                  ifelse(ncol(ctoi_mat) >= 5, round(0.2*ncol(ctoi_mat)), 1))

    if (length(random_type_samples) == 1) {
      as.vector(ctoi_mat[,random_type_samples])
    }else{
      Rfast::rowMedians(ctoi_mat[,random_type_samples], parallel = TRUE)
    }

  })
  rownames(pure_ctoi_mat) <- rownames(ref)

  # Make fractions for pure_ctoi_mat
  frac <- c(0, runif(n_samples-2, 0.01, 0.25), 0.25)
  pure_ctoi_mat.fracs <- t(t(pure_ctoi_mat) * frac)
  colnames(pure_ctoi_mat.fracs) <- make.unique(rep(ctoi, ncol(pure_ctoi_mat.fracs)))

  # Add background cell types to the pure CTOI matrix
  getFracs <- function(n_types, ctoi_frac){
    x <- runif(n_types, 0, 1)
    fracs <- x * (1-ctoi_frac) / sum(x)
    return(round(fracs, 4))
  }

  # Choose control cell types
  not_dep_cells <- celltypes[!celltypes %in% c(ctoi, dep_list[[ctoi]])] # Use only independent cell types
  pure_ct_mat.nodep <- pure_ct_mat[,colnames(pure_ct_mat) %in% not_dep_cells]
  n_types <- ifelse(length(not_dep_cells) > 7, 7, length(not_dep_cells))  # Take maximum of 7 cell types for the mixture

  # Add control cells to the corresponding columns in pure_ctoi_mat.fracs
  for (i in 1:ncol(pure_ctoi_mat.fracs)) {
    types_random_fracs <- getFracs(n_types, ctoi_frac = frac[i])     # Get control cell types fractions
    random_type_samples <- sample(1:ncol(pure_ct_mat.nodep), n_types)  # Sample n_types of the background cell types
    types_frac_vector <- rowSums(pure_ct_mat.nodep[,random_type_samples] %*% diag(types_random_fracs))  # Multiply each pure cell type by the corresponding fraction and sum by rows
    pure_ctoi_mat.fracs[,i] <- pure_ctoi_mat.fracs[,i] + types_frac_vector # Add the background cell type fractions to the fraction of CTOI
  }



  out <- list(mixture = pure_ctoi_mat.fracs,
              ctoi_fracs = frac)

  return(out)
}


# Plot signatures heatmap for a cell type
plotHeatMap <- function(type, scores_mat_tidy, signatures_collection_filtered = NULL, cor_mat){

  sig_score.df <- scores_mat_tidy %>%
    filter(signature_ct == type) %>%
    dplyr::select(signature, sample_ct, score) %>%
    pivot_wider(names_from = sample_ct, values_from = score) %>%
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






