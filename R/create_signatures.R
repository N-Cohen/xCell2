
createSignatures <- function(ref, labels, dep_list, quantiles_matrix, probs, diff_vals, cor_mat, min_genes, max_genes, extra4sim, up_genes){


  celltypes <- unique(labels[,2])
  samples <- labels[,2]
  updown <- ifelse(up_genes, "up", "down")

  all_sigs <- list()
  for (type in celltypes){

    # Get dependent cell types and remove them from the quantiles matrix
    dep_cells <- dep_list[[type]]
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, dep_cells)]

    # Should we add extra genes for similar cell-types?
    if (extra4sim) {
      # Find similar cell types
      cor_values <- cor_mat[type, ]
      cor_values <- sort(cor_values[not_dep_celltypes], decreasing = TRUE)

      # Remove cell-types with dependencies in cor_values
      to_remove <- c()
      for (i in 1:length(cor_values)) {
        if (names(cor_values[i]) %in% to_remove) {
          next
        }
        deps <- dep_list[[names(cor_values[i])]]
        to_remove <- c(to_remove, names(cor_values)[names(cor_values) %in% deps])
      }
      cor_values <- cor_values[!names(cor_values) %in% to_remove]

      cor_cells_cutoff <- quantile(cor_values, 0.9, na.rm=TRUE)
      sim_cells <- names(cor_values[which(cor_values > cor_cells_cutoff & cor_values > 0.85)])
    }



    for (diff in diff_vals) {
      for (p in 1:length(probs)) {

        # Get a list of Boolean matrices with genes that pass the quantiles criteria
        if (up_genes) {
          genes <- lapply(not_dep_celltypes, function(x){
            get(type, quantiles_matrix)[p,] > get(x, quantiles_matrix)[nrow(quantiles_matrix[[1]])-p+1,] + diff
          })
        }else{
          genes <- lapply(not_dep_celltypes, function(x){
            get(type, quantiles_matrix)[p,] + ((diff+1)*50) < get(x, quantiles_matrix)[nrow(quantiles_matrix[[1]])-p+1,]
          })
        }

        genes.mat <- matrix(unlist(genes), nrow = length(genes), byrow = TRUE,
                                 dimnames = list(not_dep_celltypes, rownames(ref)))

        sim_genes <- c()
        if (extra4sim) {
          # Find signature genes for similar cell types
          sim_genes <- c()
          n_sim <- 0

          # In case there is only one similar cell type
          if (length(sim_cells) == 1) {
            n_sim <- 1
            sim_genes <- names(which(genes.mat[rownames(genes.mat) == sim_cells,] > 0))
            if (length(sim_genes) > round(max_genes*0.25)) {
              sim_genes <- sim_genes[1:round(max_genes*0.25)]
            }
          }

          # In case there are more than one cell types
          if (length(sim_cells) > 1) {
            for (n_sim in length(sim_cells):1) {
              sim_genes <- names(which(colSums(genes.mat[rownames(genes.mat) %in% sim_cells,]) >= n_sim))
              if (length(sim_genes) > 0) {
                # Make sure sim_genes is not bigger than 1/4 of max_genes
                if (length(sim_genes) > round(max_genes*0.25)) {
                  sim_genes <- sim_genes[1:round(max_genes*0.25)]
                }
                break
              }
            }
          }
        }



        # Find up signature genes for all cell types
        n_all_celltyes <- nrow(genes.mat)
        n_min_celltypes <- ifelse(up_genes, round(nrow(genes.mat)*0.5), round(nrow(genes.mat)*0.9))
        for (n_all in n_all_celltyes:n_min_celltypes) {
          genes <- names(which(colSums(genes.mat) >= n_all))
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
            sig_name <-  paste(paste0(type, "#", updown), probs[p], diff, n_sim, n_all, sep = "_")
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


