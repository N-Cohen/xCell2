
# This function return a correlation matrix given the counts and cell types
getCellTypeCorrelation <- function(counts, samples, celltypes){

  cor_mat <- matrix(1, ncol = length(celltypes), nrow = length(celltypes), dimnames = list(celltypes, celltypes))
  n <- 2
  for (i in 1:nrow(cor_mat)) {
    for (j in n:ncol(cor_mat)) {
      celltype_i <- rownames(cor_mat)[i]
      celltype_j <- colnames(cor_mat)[j]
      cor_out <- cor(rowMedians(counts[,samples == celltype_i], na.rm = TRUE),
                     rowMedians(counts[,samples == celltype_j], na.rm = TRUE))
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
