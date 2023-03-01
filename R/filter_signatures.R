
filterSignatures <- function(pure_ct_mat, dep_list, signatures_collection, score_method, grubbs_cutoff){


  celltypes <- colnames(pure_ct_mat)

  # Score signatures
  scores_mat <- matrix(nrow = length(signatures_collection),
                       ncol = ncol(pure_ct_mat),
                       dimnames = list(names(signatures_collection), colnames(pure_ct_mat)))

  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "#"), "[", 1))

  for (type in celltypes) {

    type_signatures <- signatures_collection[type == sig_type]
    dep_cells <- dep_list[[type]]

    if (length(dep_cells) > 0) {
      types_to_use <- !colnames(scores_mat) %in% dep_cells
    }else{
      types_to_use <- rep(TRUE, ncol(scores_mat))
    }

    if (score_method == "ssgsea") {
      # Score with ssGSEA
      ssgsea_out <- GSVA::gsva(pure_ct_mat[, types_to_use], signatures_collection[type == sig_type], method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)
      scores_mat[rownames(ssgsea_out), colnames(ssgsea_out)] <- ssgsea_out

    }else if(score_method == "singscore"){
      # Score with SingScore
      sub_mix <- pure_ct_mat[,types_to_use]
      sub_mix_ranked <- singscore::rankGenes(sub_mix)
      for (i in 1:length(type_signatures)) {
        sig <- type_signatures[i]
        scores_out <- singscore::simpleScore(sub_mix_ranked, upSet = sig[[1]], centerScore = FALSE)$TotalScore
        scores_mat[which(rownames(scores_mat) == names(sig)),
                   colnames(sub_mix_ranked)] <- scores_out
      }
    }
  }

  # Make score matrix tidy
  scores_mat_tidy <- scores_mat %>%
    as_tibble(., rownames = NA) %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature, values_to = "score", names_to = "sample_ct") %>%
    separate(signature, into = "signature_ct", sep = "#", remove = FALSE, extra = "drop")%>%
    drop_na()


  # Filter by Grubb's test
  grubbs <- scores_mat_tidy %>%
    group_by(signature_ct, signature) %>%
    summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
    filter(grubbs_statistic >= quantile(grubbs_statistic, grubbs_cutoff)) %>%
    pull(signature)



  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% grubbs]

  filter_signature_out <- list("scoreMatTidy" = scores_mat_tidy, "sigCollectionFilt" = signatures_collection_filtered)

  return(filter_signature_out)
}

