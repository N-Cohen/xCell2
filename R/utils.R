library(tidyverse)

# This function return a pseudo-bulk expression matrix give a single cell RNA-Seq data.
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
