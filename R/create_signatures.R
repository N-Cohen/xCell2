
# This function return signature genes per: (# of celltypes -1) to per_cells*(# of celltypes -1)
genesOut <- function(diff_genes.mat, min_genes, max_genes, per_cells, sim_cells){

  # Find genes signature for similar cell types
  sim_genes_out <- c()
  if (!is.null(sim_cells)) {
    if (length(sim_cells) == 1) {
      sim_genes_passed <- diff_genes.mat[sim_cells,]
      sim_genes_out <- names(which(sim_genes_passed))
    }else{
      sim_cells_diff_genes.mat <- diff_genes.mat[rownames(diff_genes.mat) %in% sim_cells,]
      for (i in length(sim_cells):1) {
        sim_genes_passed <- colSums(sim_cells_diff_genes.mat) >= i
        if (sum(sim_genes_passed) > 0) {
          sim_genes_out <- names(which(sim_genes_passed))
          break
        }
      }
    }
  }


  # Find gene signature for all cell types
  n_cells <- nrow(diff_genes.mat):as.integer(per_cells*nrow(diff_genes.mat))
  genes_out <- vector(mode = "list", length = length(n_cells))
  remove_empty <- c()
  for (n in 1:length(n_cells)){
    genes_passed <- names(which(colSums(diff_genes.mat) >= n_cells[n]))

    # If not genes found - remove this n_cells value
    if (length(genes_passed) == 0) {
      remove_empty <- c(remove_empty, n)
      next
    }

    # If genes_passed >= max_genes - take the first max_genes genes
    if (length(genes_passed) >= max_genes) {
      genes_out[n][[1]] <- genes_passed[1:max_genes]
      next
    }

    # Add genes of similar cell types
    all_genes_passed <- unique(c(genes_passed, sim_genes_out))

    # If all_genes_passed >= max_genes - take the first max_genes genes
    if (length(all_genes_passed) >= max_genes) {
      genes_out[n][[1]] <- all_genes_passed[1:max_genes]
      next
    }

    # If not enough genes - remove this n_cells value
    if (length(all_genes_passed) < min_genes) {
      remove_empty <- c(remove_empty, n)
      next
    }

    genes_out[n][[1]] <- all_genes_passed

  }

  # Remove # of cell type without any genes found
  if (!is.null(remove_empty)) {
    genes_out <- genes_out[-remove_empty]
    names(genes_out) <- as.character(n_cells)[-remove_empty]
  }else{
    names(genes_out) <- as.character(n_cells)
  }

  return(genes_out)
}

# This function return signature genes per # of celltypes per quantiles
quantilesGenesOut <- function(quantiles_matrix, probs, type, diff, min_genes, max_genes, per_cells, sim_cells){

  quantiles_genes_out <- vector(mode = "list", length = as.integer(length(probs)/2 + 0.5))
  for (j in 1:as.integer(length(probs)/2 + 0.5)) {
    names(quantiles_genes_out)[j] <- paste0(round(probs[j], 2)*100, "%,",
                                            round(probs[length(probs)-(j-1)], 2)*100, "%")

    # Get diff_genes for this diff and quantile values
    type_id <- which(names(quantiles_matrix) == type)
    diff_genes <- lapply(quantiles_matrix[-type_id],
                         function(x) quantiles_matrix[[type_id]][j,] > x[length(probs)-(j-1),]+diff)
    diff_genes.mat <- matrix(unlist(diff_genes), nrow = length(diff_genes), byrow = TRUE)
    colnames(diff_genes.mat) <- colnames(quantiles_matrix[[type_id]])
    rownames(diff_genes.mat) <- names(quantiles_matrix)[-type_id]

    # Select genes
    genes_out <- genesOut(diff_genes.mat, min_genes, max_genes, per_cells, sim_cells)
    if (length(genes_out) == 0) {
      next
    }
    quantiles_genes_out[j][[1]] <- genes_out
  }

  # If all quantiles values are empty - return NULL
  quantiles_genes_out <- Filter(Negate(is.null), quantiles_genes_out)
  if (length(quantiles_genes_out) == 0) {
    return(NULL)
  }
  return(quantiles_genes_out)
}


diffQuantilesGenesOut <- function(diff_vals, quantiles_matrix, probs, min_genes, max_genes, per_cells, type, sim_cells){

  diff_quantiles_genes_out <- vector(mode = "list", length = length(diff_vals))
  for (i in 1:length(diff_vals)) {
    quantilesGenes_out <- quantilesGenesOut(quantiles_matrix, probs, type, diff_vals[i], min_genes, max_genes,
                                            per_cells, sim_cells)

    if (length(quantilesGenes_out) == 0) {
      next
    }

    diff_quantiles_genes_out[i][[1]] <- quantilesGenes_out
    names(diff_quantiles_genes_out)[i] <- as.character(diff_vals[i])
  }

  diff_quantiles_genes_out <- Filter(Negate(is.null), diff_quantiles_genes_out)

  if (length(diff_quantiles_genes_out) == 0) {
    return(NULL)
  }
  return(diff_quantiles_genes_out)
}


# Load this for debugging (Remove):
diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
probs = c(.1, .25, .33333333, .5, .6666666, .75, .9)
min_genes = 7; max_genes = 200; per_cells = .95; cor_cells_cutoff = .92
type_id = 30; type = "CL:0000057"
quantiles_matrix <- readRDS("quantiles_matrix.RData")



createSignatures <- function(counts, samples, celltypes, dep_list = NULL, celltype_cor_mat = NULL,
                             diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                             probs = c(.1, .25, .33333333, .5, .6666666, .75, .9),
                             min_genes = 7, max_genes = 200, per_cells = .95, cor_cells_cutoff = .92){


  # (1) Calculate quantiles for each cell type
  message("Calculating quantiles...")

  quantiles_matrix <- lapply(celltypes, function(type){
    # If there is one sample for this cell type - duplicate the sample to make a data frame
    if (sum(samples == type, na.rm = TRUE) == 1) {
      type.df <- cbind(counts[,samples==type & !is.na(samples)], counts[,samples==type & !is.na(samples)])
    }else{
      type.df <- counts[,samples==type & !is.na(samples)]
    }
    quantiles_matrix <- apply(type.df, 1, function(x) quantile(x, probs, na.rm=TRUE))
  })
  names(quantiles_matrix) <- celltypes


  # (2) Create signature for each cell type
  message("Finding signatures...")

  signatures_out <- vector(mode = "list", length = length(celltypes))
  for (type_id in 1:length(celltypes)) {
    type <- celltypes[type_id]
    names(signatures_out)[type_id] <- type
    message(type)

    # Find dependent cell types
    dep_cells <- NULL
    if (!is.null(dep_list)) {
      dep_cells <- dep_list[[type]]
      if (length(dep_cells) > 0) {
        quantiles_matrix_no_deps <- quantiles_matrix[!names(quantiles_matrix) %in% dep_cells]
      }else{
        quantiles_matrix_no_deps <- quantiles_matrix
      }
    }else{
      quantiles_matrix_no_deps <- quantiles_matrix
    }


    # Find similar cell types
    sim_cells <- NULL
    if (!is.null(celltype_cor_mat)){
      sim_cells <- names(which(celltype_cor_mat[type, ] >= cor_cells_cutoff))
      sim_cells <- sim_cells[!sim_cells %in% c(type, dep_cells)]
      if (length(sim_cells) == 0) {
        sim_cells <- NULL
      }
    }

    diffQuantilesGenes_out <- diffQuantilesGenesOut(diff_vals, quantiles_matrix_no_deps,
                                                    probs, min_genes, max_genes,
                                                    per_cells, type, sim_cells)

    if (length(diffQuantilesGenes_out) == 0) {
      next
    }

    signatures_out[type_id][[1]] <- diffQuantilesGenes_out
  }

  if (length(signatures_out) == 0) {
    warning("No signatures found for reference")
  }


  return(signatures_out)
}



# Testing -----------------------------

# hpca <- readRDS("~/Documents/xCell2.0/HPCA.RData")
# ref_list <- list("HPCA" = hpca)

# hpca.out <- createSignatures(ref_list, RNA_seq = FALSE, cor_cells_cutoff = 0.92, OBOfile = "~/Documents/xCell2.0/cl.obo")
# saveRDS(hpca.out, file="~/Documents/xCell2.0/HPCA_signatures.RData")





# library(singscore)
# library(tidyverse)
# library(pheatmap)
#
# HPCAsig <- readRDS("~/Documents/xCell2.0/HPCA_signatures.RData")
# # HPCAsig <- hpca.out
#
# # Rank reference
# hpca <- readRDS("~/Documents/xCell2.0/HPCA.RData")
# hpca.ranked <- rankGenes(hpca)
#
# # Make heatmap for all signatures of CL:0000840
# example_signature <- HPCAsig$HPCA$`CL:0000451`$`0`$`10%,90%`$`62`
#
# score.df <- simpleScore(hpca.ranked, upSet = example_signature)
# df <- score.df %>%
#   select(TotalScore) %>%
#   mutate(celltype = hpca$label.ont) %>%
#   drop_na() %>%
#   group_by(celltype) %>%
#   summarise(sig1 = median(TotalScore)) %>%
#   arrange(-sig1)
#
# celltype.list <- HPCAsig$HPCA$`CL:0000840`
#
#
# for (d in 1:length(celltype.list)) {
#   diff <- names(celltype.list)[d]
#   for (q in 1:length(celltype.list[[d]])) {
#     quantile <- names(celltype.list[[d]])[q]
#     for (n in 1:length(celltype.list[[d]][[q]])) {
#       ncells <- names(celltype.list[[d]][[q]])[n]
#
#       sig <- celltype.list[[d]][[q]][[n]]
#       sig_name <- paste0(diff, ";", quantile, ";", ncells, ";", length(sig))
#       print(sig_name)
#
#
#       score.df <- simpleScore(hpca.ranked, upSet = sig)
#
#       df <- score.df %>%
#         select(TotalScore) %>%
#         mutate(celltype = hpca$label.ont) %>%
#         drop_na() %>%
#         group_by(celltype) %>%
#         summarise(medianScore = median(TotalScore)) %>%
#         rename(!!sig_name := medianScore) %>%
#         full_join(df)
#     }
#   }
# }
#
# mat <- df %>%
#   select(-sig1) %>%
#   tibble::column_to_rownames(var = "celltype") %>%
#   as.matrix()
#
#
# pheatmap(t(mat))
#
#
# good_signature <- HPCAsig$HPCA$`CL:0000840`$`4`$`33%,67%`$`60`
# score.df <- simpleScore(hpca.ranked, upSet = good_signature)
#
# score.df %>%
#   mutate(ont = hpca$label.ont, fine = hpca$label.fine, main = hpca$label.main) %>%
#   select(-TotalDispersion) %>%
#   ggplot(.,aes(x=reorder(ont, -TotalScore), y=TotalScore)) +
#   geom_boxplot() +
#   geom_point() +
#   labs(y="Signature Score", x="") +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
