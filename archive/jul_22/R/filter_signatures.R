

filterSignatures <- function(signatures_collection, scores_mat_pure_tidy, take_top_per){

  # Fixing bug (part 1): if a cell type do not have > 2 scores with other cell types (because he is the mixtures control)
  lonely_celltype <- scores_mat_pure_tidy %>%
    drop_na() %>%
    group_by(signature_ct) %>%
    summarise(n = length(unique(mixture_ct))) %>%
    filter(n <= 2) %>%
    pull(signature_ct)


  # Rank signature by Grubbs' test
  signatures_filtered <- scores_mat_pure_tidy %>%
    filter(!signature_ct %in% lonely_celltype) %>%
    drop_na() %>%
    # Remove signatures which the ctoi score is not the max_score
    group_by(signature_ct, signature) %>%
    filter(any(signature_ct == mixture_ct & score == max(score))) %>%
    # Rank signatures
    summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$p.value) %>%
    mutate(grubbs_rank = percent_rank(dplyr::desc(grubbs_pvalue))*100) %>%
    arrange(-grubbs_rank, .by_group = TRUE) %>%
    # Filter signatures
    filter(grubbs_rank >= quantile(grubbs_rank, 1-take_top_per, na.rm = TRUE))


  # Maximum 20 signature per cell type
  ct_to_filt <- names(table(signatures_filtered$signature_ct)[table(signatures_filtered$signature_ct) > 20])
  ct_to_keep <- signatures_filtered %>%
    filter(!signature_ct %in% ct_to_filt)

  signatures_filtered <- signatures_filtered %>%
    group_by(signature_ct) %>%
    filter(signature_ct %in% ct_to_filt) %>%
    filter(row_number() %in% 1:20) %>%
    rbind(ct_to_keep)


  # Fixing bug (part 2)
  if (length(lonely_celltype) > 0) {
    lonely_celltype_signatures_filtered <- scores_mat_pure_tidy %>%
      filter(signature_ct %in% lonely_celltype) %>%
      drop_na() %>%
      group_by(signature_ct, signature) %>%
      filter(any(signature_ct == mixture_ct & score == max(score))) %>%
      mutate(order = ifelse(signature_ct == mixture_ct, 2, 1)) %>%
      summarise(delta = score - lag(score, order_by = order)) %>%
      drop_na() %>%
      # Filter signatures
      group_by(signature_ct) %>%
      filter(delta >= quantile(delta, 1-take_top_per, na.rm = TRUE)) %>%
      dplyr::select(-delta) %>%
      mutate(signature_ct = lonely_celltype, grubbs_pvalue = NA, grubbs_rank = NA) %>%
      dplyr::select(signature_ct, signature, grubbs_pvalue, grubbs_rank)

    signatures_filtered <- rbind(signatures_filtered, lonely_celltype_signatures_filtered)
  }

  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% signatures_filtered$signature]


  return(signatures_collection_filtered)
}

#
