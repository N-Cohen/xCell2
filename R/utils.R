
# This function return a correlation matrix given the counts and cell types
getCellTypeCorrelation <- function(counts, samples, celltypes){

  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  n <- 2
  for (i in 1:nrow(cor_mat)) {
    for (j in n:ncol(cor_mat)) {
      celltype_i <- rownames(cor_mat)[i]
      celltype_j <- colnames(cor_mat)[j]

      if (sum(samples == celltype_i) == 1) {
        median_expression_i <- as.numeric(counts[,samples == celltype_i])
      }else{
        median_expression_i <- rowMedians(counts[,samples == celltype_i], na.rm = TRUE)
      }

      if (sum(samples == celltype_j) == 1) {
        median_expression_j <- as.numeric(counts[,samples == celltype_j])
      }else{
        median_expression_j <- rowMedians(counts[,samples == celltype_j], na.rm = TRUE)
      }

      cor_out <- cor(median_expression_i, median_expression_j)
      cor_mat[i, j] <- cor_out
      cor_mat[j, i] <- cor_out
    }
    n <- n + 1
    if (n > ncol(cor_mat)) {
      break
    }
  }

  return(cor_mat)
}


# This function return a vector of cell type dependencies give a list of cell types and one cell type on interest (type)
getDependencies <- function(OBOfile, celltypes){

  cl <- suppressWarnings(ontologyIndex::get_ontology(OBOfile))
  celltype2dep <- vector(mode = "list", length = length(celltypes))
  names(celltype2dep) <- celltypes

  for (type in celltypes) {
    descendants <- ontologyIndex::get_descendants(cl, roots = type, exclude_roots = T)
    ancestors <- ontologyIndex::get_ancestors(cl, terms = type)
    ancestors <- ancestors[ancestors != type]
    dep_cells <- celltypes[celltypes %in% c(descendants, ancestors)]
    celltype2dep[[type]] <- dep_cells
  }

  return(celltype2dep)
}

#  TODO: add options for synthetic control
createInSilicoMixture <- function(counts, celltypes, celltype_cor_mat,
                                  synthetic_control = FALSE, fractions = seq(.01,.25, .01)){


  fractions_mat_list <- lapply(celltypes, function(type) {

    # Get control vector
    control <- names(which.min(celltype_cor_mat[type,]))
    control_samples <- colnames(ref)[which(ref$label.ont == control)]
    if (length(control_samples) == 1) {
      control_vec <- as.vector(counts[,colnames(counts) == control_samples])
    }else{
      control_vec <- as.vector(apply(counts[,colnames(counts) %in% control_samples], 1, median))
    }

    # Get CTOI vector
    type_samples <- colnames(ref)[which(ref$label.ont == type)]
    if (length(type_samples) == 1) {
      type_vec <- as.vector(counts[,colnames(counts) %in% type_samples])
    }else{
      type_vec <- as.vector(apply(counts[,colnames(counts) %in% type_samples], 1, median))
    }

    # Calculate fractions
    fractions <- c(1, fractions) # For pure CTOI scoring
    fractions_mat <- sapply(fractions, function(f) {
      type_vec * f + control_vec*(1-f)
    })
    colnames(fractions_mat) <- paste(fractions, control, sep = "_")
    fractions_mat
  })

  names(fractions_mat_list) <- celltypes
  ref_insilico_mat <- do.call(cbind.data.frame, fractions_mat_list)
  colnames(ref_insilico_mat) <- sub("\\.", "_", colnames(ref_insilico_mat))
  rownames(ref_insilico_mat) <- rownames(ref)
  return(ref_insilico_mat)
}





# Helper function (for development) ------


# Plot signatures heatmap for a cell type
plotHeatMap <- function(ctoi, scores_mat, dep_list, ref, signatures_ranked, ranks_weights = c(1, 1, 1, 1), take_top_per = NULL, label = "ont"){

  chnageLabel <- function(ref, ont, label = "fine"){
    if (label == "main") {
      out <- unique(ref$label.main[ref$label.ont == ont])
      return(out[!is.na(out)][1])
    }else if(label == "fine"){
      out <- unique(ref$label.fine[ref$label.ont == ont])
      return(out[!is.na(out)][1])
    }
  }

  # Take only scores with 100% CTOI mixture
  scores_mat_sub <- scores_mat[startsWith(rownames(scores_mat), ctoi), grep("*_1_", colnames(scores_mat))]
  colnames(scores_mat_sub) <- unlist(lapply(strsplit(colnames(scores_mat_sub), "_"), "[", 1))

  # Remove dep celltypes
  types2use <- names(sort(celltype_cor_mat[ctoi,], decreasing = T))
  types2use <- types2use[!types2use %in% dep_list[[ctoi]]]

  # Sort signatures by rank
  sigs_sorted <- signatures_ranked %>%
    filter(signature_ct == ctoi) %>%
    mutate(CT_rank = coalesce(CT_rank, CTOI_rank), diff_rank = coalesce(diff_rank, CTOI_rank)) %>%
    mutate(final_rank = ranks_weights[1]*CTOI_rank + ranks_weights[2]*CT_rank + ranks_weights[3]*diff_rank + ranks_weights[4]*grubbs_rank) %>%
    group_by(signature_ct) %>%
    arrange(-final_rank, .by_group = TRUE)

  # Filter signatures
  if (!is.null(take_top_per)) {
    sigs_sorted <- sigs_sorted %>%
      filter(final_rank > quantile(final_rank, 1-take_top_per))
  }

  scores_mat_sub <- scores_mat_sub[sigs_sorted$signature, types2use]

  # Change celltype labels
  if (label == "main") {
    colnames(scores_mat_sub) <- unlist(lapply(colnames(scores_mat_sub), function(x) chnageLabel(ref, x, label = "main")))
  }else if(label == "fine"){
    colnames(scores_mat_sub) <- unlist(lapply(colnames(scores_mat_sub), function(x) chnageLabel(ref, x)))
  }

  hm <- pheatmap::pheatmap(scores_mat_sub, cluster_rows=F, cluster_cols=F, scale = "row", col=brewer.pal(11, "RdBu")) # Change the yellow color
  return(hm)
}


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

