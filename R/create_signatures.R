
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, cor_mat, diff_vals, min_genes, max_genes){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  all_sigs <- list()
  for (type in celltypes){
    # Get dependent cell types and remove them from the quantiles matrix.
    dep_cells <- dep_list[[type]]
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, dep_cells)]

    # Find similar cell types
    cor_values <- cor_mat[type, ]
    cor_values <- sort(cor_values[!names(cor_values) %in% c(type, dep_cells)], decreasing = TRUE)

    # Remove cell-types with dependencies in cor_values (new!)
    to_remove <- c()
    for (i in 1:length(cor_values)) {
      if (names(cor_values[i]) %in% to_remove) {
        next
      }
      deps <- dep_list[[names(cor_values[i])]]
      to_remove <- c(to_remove, names(cor_values)[names(cor_values) %in% deps])
    }
    cor_values <- cor_values[!names(cor_values) %in% to_remove]

    # TODO: Should we add extra genes for similar cell-types or not? What are similar cell-types?
    cor_cells_cutoff <- quantile(cor_values, 0.9, na.rm=TRUE)
    sim_cells <- names(cor_values[which(cor_values > cor_cells_cutoff & cor_values > 0.85)])
    # sim_cells <- names(cor_values[which(cor_values > 0.5)])
    # sim_cells <- c()

    for (diff in diff_vals) {
      for (p in 1:length(probs)) {

        # Get a list of Boolean matrices with genes that pass the quantiles criteria
        diff_genes <- lapply(not_dep_celltypes, function(x){
          get(type, quantiles_matrix)[p,] > get(x, quantiles_matrix)[nrow(quantiles_matrix[[1]])-p+1,] + diff
        })

        diff_genes.mat <- matrix(unlist(diff_genes), nrow = length(diff_genes), byrow = TRUE,
                                 dimnames = list(not_dep_celltypes, rownames(ref)))


        # Find signature genes for similar cell types
        sim_genes <- c()
        n_sim <- 0

        # In case there is only one similar cell type
        if (length(sim_cells) == 1) {
          n_sim <- 1
          sim_genes <- names(which(diff_genes.mat[rownames(diff_genes.mat) == sim_cells,] > 0))
          if (length(sim_genes) > round(max_genes*0.25)) {
            sim_genes <- sim_genes[1:round(max_genes*0.25)]
          }
        }

        # In case there are more than one cell types
        if (length(sim_cells) > 1) {
          for (n_sim in length(sim_cells):1) {
            sim_genes <- names(which(colSums(diff_genes.mat[rownames(diff_genes.mat) %in% sim_cells,]) >= n_sim))
            if (length(sim_genes) > 0) {
              # Make sure sim_genes is not bigger than 1/4 of max_genes
              if (length(sim_genes) > round(max_genes*0.25)) {
                sim_genes <- sim_genes[1:round(max_genes*0.25)]
              }
              break
            }
          }
        }


        # Find signature genes for all cell types
        for (n_all in nrow(diff_genes.mat):round(nrow(diff_genes.mat)*0.5)) {
          genes <- names(which(colSums(diff_genes.mat) >= n_all))
          # Merge with sim_genes
          genes <- unique(c(sim_genes, genes))
          # If check there are enough genes
          if (length(genes) < min_genes) {
            next
          }
          # Check if there are too many genes
          if (length(genes) > max_genes) {
            genes <- genes[1:max_genes]
          }
          # Save signature
          if (length(genes) > 0) {
            sig_name <-  paste(paste0(type, "#"), probs[p], diff, n_sim, n_all, sep = "_")
            all_sigs[[sig_name]] <- GSEABase::GeneSet(genes, setName = sig_name)
          }
        }
      }
    }
  }

  if (length(all_sigs) == 0) {
    warning("No signatures found for reference!")
  }

  # Make GeneSetCollection object
  signatures_collection <- GSEABase::GeneSetCollection(all_sigs)

  return(signatures_collection)
}

