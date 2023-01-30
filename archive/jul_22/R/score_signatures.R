
scoreSignatures <- function(signatures_collection, inSilico_mat, dep_list, score_method){

  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "_"), "[", 1))
  controls <- unlist(lapply(strsplit(colnames(inSilico_mat), "_"), "[", 3))
  mixture_type <- unlist(lapply(strsplit(colnames(inSilico_mat), "_"), "[", 1))

  # Make an empty score_mat matrix
  scores_mat <- matrix(nrow = length(signatures_collection), ncol = ncol(inSilico_mat), dimnames = list(names(signatures_collection), colnames(inSilico_mat)))

  for (type in unique(sig_type)) {

    type_signatures <- signatures_collection[type == sig_type]
    dep_cells <- dep_list[[type]]

    # Use only mixtures in which: (1) The mixture's type not in dep_cells
    #                             (2) The mixture's control not in dep_cells
    #                             (2) The mixture's control not the CTOI
    mixtures_to_use <- (mixture_type %in% dep_cells) + (controls %in% dep_cells) + (controls == type) == 0

    # TODO: WARNING!! need to fix this for ABIS "B cells" type
    # sum(mixtures_to_use) == 1 because both the dep cells and controls limit the available in silico mixture
    if (sum(mixtures_to_use) == 1) {
      mixtures_to_use <- (mixture_type %in% dep_cells) + (controls == type) == 0
    }


    if (score_method == "ssgsea") {
      # Score with ssGSEA
      ssgsea_out <- GSVA::gsva(inSilico_mat[, mixtures_to_use], signatures_collection[type == sig_type], method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)
      scores_mat[rownames(ssgsea_out), colnames(ssgsea_out)] <- ssgsea_out

    }else if(score_method == "singscore"){
      # Score with SingScore
      sub_mix <- inSilico_mat[,mixtures_to_use]
      ctoi_ranked_inSilico_mat <- singscore::rankGenes(sub_mix)
      for (i in 1:length(type_signatures)) {
        sig <- type_signatures[i]
        scores_out <- singscore::simpleScore(ctoi_ranked_inSilico_mat, upSet = sig[[1]], centerScore = FALSE)$TotalScore
        scores_mat[which(rownames(scores_mat) == names(sig)),
                   colnames(ctoi_ranked_inSilico_mat)] <- scores_out
      }
    }
  }


  return(scores_mat)
}


# Testing ---------------------------------------
