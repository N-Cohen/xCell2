#
# # This function return signature genes per: (# of celltypes -1) to per_cells*(# of celltypes -1)
# genesOut <- function(diff_genes.mat, min_genes, max_genes, per_cells, sim_cells){
#
#   # Find genes signature for similar cell types
#   sim_genes_out <- c()
#   if (!is.null(sim_cells)) {
#     if (length(sim_cells) == 1) {
#       sim_genes_passed <- diff_genes.mat[sim_cells,]
#       sim_genes_out <- names(which(sim_genes_passed))
#     }else{
#       sim_cells_diff_genes.mat <- diff_genes.mat[rownames(diff_genes.mat) %in% sim_cells,]
#       for (i in length(sim_cells):1) {
#         sim_genes_passed <- colSums(sim_cells_diff_genes.mat) >= i
#         if (sum(sim_genes_passed) > 0) {
#           sim_genes_out <- names(which(sim_genes_passed))
#           break
#         }
#       }
#     }
#   }
#
#
#   # Find gene signature for all cell types
#   n_cells <- nrow(diff_genes.mat):as.integer(per_cells*nrow(diff_genes.mat))
#   genes_out <- vector(mode = "list", length = length(n_cells))
#   remove_empty <- c()
#   for (n in 1:length(n_cells)){
#     genes_passed <- names(which(colSums(diff_genes.mat) >= n_cells[n]))
#
#     # If not genes found - remove this n_cells value
#     if (length(genes_passed) == 0) {
#       remove_empty <- c(remove_empty, n)
#       next
#     }
#
#     # If genes_passed >= max_genes - take the first max_genes genes
#     if (length(genes_passed) >= max_genes) {
#       genes_out[n][[1]] <- genes_passed[1:max_genes]
#       next
#     }
#
#     # Add genes of similar cell types
#     all_genes_passed <- unique(c(genes_passed, sim_genes_out))
#
#     # If all_genes_passed >= max_genes - take the first max_genes genes
#     if (length(all_genes_passed) >= max_genes) {
#       genes_out[n][[1]] <- all_genes_passed[1:max_genes]
#       next
#     }
#
#     # If not enough genes - remove this n_cells value
#     if (length(all_genes_passed) < min_genes) {
#       remove_empty <- c(remove_empty, n)
#       next
#     }
#
#     genes_out[n][[1]] <- all_genes_passed
#
#   }
#
#   # Remove # of cell type without any genes found
#   if (!is.null(remove_empty)) {
#     genes_out <- genes_out[-remove_empty]
#     names(genes_out) <- as.character(n_cells)[-remove_empty]
#   }else{
#     names(genes_out) <- as.character(n_cells)
#   }
#
#   return(genes_out)
# }
#
# # This function return signature genes per # of celltypes per quantiles
# quantilesGenesOut <- function(quantiles_matrix, probs, type, diff, min_genes, max_genes, per_cells, sim_cells){
#
#   quantiles_genes_out <- vector(mode = "list", length = as.integer(length(probs)/2 + 0.5))
#   for (j in 1:as.integer(length(probs)/2 + 0.5)) {
#     names(quantiles_genes_out)[j] <- paste0(round(probs[j], 2)*100, "%,",
#                                             round(probs[length(probs)-(j-1)], 2)*100, "%")
#
#     # Get diff_genes for this diff and quantile values
#     type_id <- which(names(quantiles_matrix) == type)
#     diff_genes <- lapply(quantiles_matrix[-type_id],
#                          function(x) quantiles_matrix[[type_id]][j,] > x[length(probs)-(j-1),]+diff)
#     diff_genes.mat <- matrix(unlist(diff_genes), nrow = length(diff_genes), byrow = TRUE)
#     colnames(diff_genes.mat) <- colnames(quantiles_matrix[[type_id]])
#     rownames(diff_genes.mat) <- names(quantiles_matrix)[-type_id]
#
#     # Select genes
#     genes_out <- genesOut(diff_genes.mat, min_genes, max_genes, per_cells, sim_cells)
#     if (length(genes_out) == 0) {
#       next
#     }
#     quantiles_genes_out[j][[1]] <- genes_out
#   }
#
#   # If all quantiles values are empty - return NULL
#   quantiles_genes_out <- Filter(Negate(is.null), quantiles_genes_out)
#   if (length(quantiles_genes_out) == 0) {
#     return(NULL)
#   }
#   return(quantiles_genes_out)
# }
#
#
# diffQuantilesGenesOut <- function(diff_vals, quantiles_matrix, probs, min_genes, max_genes, per_cells, type, sim_cells){
#
#   diff_quantiles_genes_out <- vector(mode = "list", length = length(diff_vals))
#
#   for (i in 1:length(diff_vals)) {
#
#     quantilesGenes_out <- quantilesGenesOut(quantiles_matrix, probs, type, diff_vals[i], min_genes, max_genes,
#                                             per_cells, sim_cells)
#
#     if (length(quantilesGenes_out) == 0) {
#       next
#     }
#
#     diff_quantiles_genes_out[i][[1]] <- quantilesGenes_out
#     names(diff_quantiles_genes_out)[i] <- as.character(diff_vals[i])
#   }
#
#   diff_quantiles_genes_out <- Filter(Negate(is.null), diff_quantiles_genes_out)
#
#   if (length(diff_quantiles_genes_out) == 0) {
#     return(NULL)
#   }
#   return(diff_quantiles_genes_out)
# }


# Load this for debugging (Remove):
#type = celltypes[4]; diff = 1; p = 1



createSignatures <- function(counts, labels, quantiles_matrix, probs, cor_mat, diff_vals, per_cells, genes_type, min_genes, max_genes){

  celltypes <- unique(labels[,2])
  samples <- labels[,2]

  all_sigs <- list()
  for (type in celltypes){
    # Get dependent cell types and remove them from the quantiles matrix
    dep_cells <- dep_list[[type]]
    not_dep_celltypes <- celltypes[!celltypes %in% c(type, dep_cells)]

    # Find similar cell types
    sim_cells <- names(which(cor_mat[type, ] >= cor_cells_cutoff))
    sim_cells <- sim_cells[!sim_cells %in% c(type, dep_cells)]


    for (diff in diff_vals) {
      for (p in 1:length(probs)) {

        # Get a list of Boolean matrices with genes that pass the quantiles criteria
        diff_genes <- lapply(not_dep_celltypes, function(x){
          get(type, quantiles_matrix)[p,] > get(x, quantiles_matrix)[nrow(quantiles_matrix[[1]])-p+1,] + diff
        })

        diff_genes.mat <- matrix(unlist(diff_genes), nrow = length(diff_genes), byrow = TRUE,
                                 dimnames = list(not_dep_celltypes, rownames(counts)))


        # Find signature genes for similar cell types
        sim_genes <- c()
        n_sim <- 0
        if (length(sim_cells) > 1) { # Minimum number of sim_cells is 2
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
        for (n_all in nrow(diff_genes.mat):round(nrow(diff_genes.mat)*per_cells)) {
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
          sig_name <-  paste(type, probs[p], diff, n_sim, n_all, sep = "_")
          all_sigs[[sig_name]] <- GeneSet(genes, setName = sig_name, geneIdType = genes_type)

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






#   # --- old ----------------------------------------------
#   signatures_out <- vector(mode = "list", length = length(celltypes))
#   for (type_id in 1:length(celltypes)) {
#     type <- celltypes[type_id]
#     names(signatures_out)[type_id] <- type
#     message(type)
#
#     # Find dependent cell types
#     dep_cells <- NULL
#     if (!is.null(dep_list)) {
#       dep_cells <- dep_list[[type]]
#       if (length(dep_cells) > 0) {
#         quantiles_matrix_no_deps <- quantiles_matrix[!names(quantiles_matrix) %in% dep_cells]
#       }else{
#         quantiles_matrix_no_deps <- quantiles_matrix
#       }
#     }else{
#       quantiles_matrix_no_deps <- quantiles_matrix
#     }
#
#
#     # Find similar cell types
#     sim_cells <- NULL
#     if (!is.null(celltype_cor_mat)){
#       sim_cells <- names(which(celltype_cor_mat[type, ] >= cor_cells_cutoff))
#       sim_cells <- sim_cells[!sim_cells %in% c(type, dep_cells)]
#       if (length(sim_cells) == 0) {
#         sim_cells <- NULL
#       }
#     }
#
#     diffQuantilesGenes_out <- diffQuantilesGenesOut(diff_vals, quantiles_matrix_no_deps,
#                                                     probs, min_genes, max_genes,
#                                                     per_cells, type, sim_cells)
#
#     if (length(diffQuantilesGenes_out) == 0) {
#       next
#     }
#
#     signatures_out[type_id][[1]] <- diffQuantilesGenes_out
#   }
#
#   if (length(signatures_out) == 0) {
#     warning("No signatures found for reference")
#   }
#
#
#   return(signatures_out)
# }



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
