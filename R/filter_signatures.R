
filterSignatures <- function(pure_ct_mat, dep_list, signatures_collection, score_method, take_top_per, max_sigs){


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
    separate(signature, into = "signature_ct", sep = "#", remove = FALSE, extra = "drop")

  # Remove signatures which the ctoi score is not the max_score
  signatures_filtered <- scores_mat_tidy %>%
    drop_na() %>%
    group_by(signature_ct, signature) %>%
    filter(any(signature_ct == sample_ct & score == max(score)))

  # Take all signatures of cell types that did not pass the filtering above
  cts_did_not_pass <- celltypes[!celltypes %in% unique(signatures_filtered$signature_ct)]
  warning(paste("Poor signature scores are found in some cell types:", paste(cts_did_not_pass, collapse = ", ")))
  for (ct in cts_did_not_pass) {
    ct_scores <- scores_mat_tidy %>%
      filter(signature_ct == ct & signature_ct == sample_ct) %>%
      drop_na() %>%
      dplyr::rename(score_ct = score)

    signatures_filtered <- rbind(signatures_filtered, scores_mat_tidy %>%
                                   filter(signature_ct == ct) %>%
                                   drop_na() %>%
                                   left_join(dplyr::select(ct_scores, signature, score_ct), by = "signature") %>%
                                   rowwise() %>%
                                   filter(score <= score_ct) %>%
                                   dplyr::select(-score_ct))
  }


  # Filter signature by Grubbs' test
  signatures_filtered.grubbs <- signatures_filtered %>%
    # Rank signatures
    summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$p.value) %>%
    mutate(grubbs_rank = percent_rank(dplyr::desc(grubbs_pvalue))*100) %>%
    arrange(-grubbs_rank, .by_group = TRUE) %>%
    # Filter signatures
    filter(if (n() >= 10) grubbs_rank >= quantile(grubbs_rank, 1-take_top_per, na.rm = TRUE) else grubbs_rank >= 0) %>% # Minimum ten signatures for each cell type
    pull(signature)


  # # Top signatures by CTOI delta score
  # signatures_filtered %>%
  #   filter(signature %in% signatures_filtered.grubbs) %>%
  #   filter(signature_ct == sample_ct) %>%
  #   filter(fraction %in% c(min(mixture_fractions[mixture_fractions != 0]), max(mixture_fractions))) %>%
  #   group_by(signature, sample_ct) %>%
  #   summarise(delta = score - lag(score)) %>%
  #   drop_na() %>%
  #   group_by(sample_ct) %>%
  #   arrange(desc(delta)) %>%
  #   top_n(3) %>%
  #   separate(signature, into = "signature_ct", sep = "#", remove = FALSE, extra = "drop")


  # # Maximum max_sigs signature per cell type
  #  ct_to_filt <- names(table(signatures_filtered$signature_ct)[table(signatures_filtered$signature_ct) > max_sigs])
  #  ct_to_keep <- signatures_filtered %>%
  #    filter(!signature_ct %in% ct_to_filt)
  #
  #  signatures_filtered <- signatures_filtered %>%
  #    group_by(signature_ct) %>%
  #    filter(signature_ct %in% ct_to_filt) %>%
  #    filter(row_number() %in% 1:max_sigs) %>%
  #    rbind(ct_to_keep)


  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% signatures_filtered.grubbs]

  filter_signature_out <- list("scoreMatTidy" = scores_mat_tidy, "sigCollectionFilt" = signatures_collection_filtered)

  return(filter_signature_out)
}

