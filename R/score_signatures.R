
scoreSignatures <- function(signatures_list, inSilico_mat, dep_list){
  celltypes <- names(signatures_list)
  all_signatures <- lapply(rapply(signatures_list, enquote, how="unlist"), eval)
  controls <- unlist(lapply(strsplit(colnames(inSilico_mat), "_"), "[", 3))
  mixture_type <- unlist(lapply(strsplit(colnames(inSilico_mat), "_"), "[", 1))

  scores_mat <- matrix(nrow = length(all_signatures), ncol = ncol(inSilico_mat), dimnames = list(names(all_signatures), colnames(inSilico_mat)))

  for (type in celltypes) {
    all_celltype_signatures <- all_signatures[grepl(paste0(type, "."), names(all_signatures))]
    dep_cells <- dep_list[[type]]

    # Use only mixtures in which: (1) control not in dep_cells (2) contol not the same as the CTOI (3) mixture_type not in dep_cells
    mixtures_to_use <- (controls %in% dep_cells) + (controls == type) + (mixture_type %in% dep_cells) == 0
    ctoi_ranked_inSilico_mat <- singscore::rankGenes(inSilico_mat[,mixtures_to_use])

    for (i in 1:length(all_celltype_signatures)) {
      sig <- all_celltype_signatures[i]
      scores_out <- singscore::simpleScore(ctoi_ranked_inSilico_mat, upSet = sig[[1]])$TotalScore
      scores_mat[names(sig), colnames(ctoi_ranked_inSilico_mat)] <- scores_out
    }
  }

  return(scores_mat)
}


