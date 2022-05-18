library(tidyverse)


rankSignatures <- function(scores_mat){

  scores_mat_tidy <- scores_mat %>%
    as_tibble(., rownames = NA) %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature,
                 names_to = c("mixture_ct", "mixture_fraction", "mixture_control"),
                 names_sep = "_",
                 values_to = "score") %>%
    filter(mixture_fraction == 1) %>%
    drop_na() %>%
    # TODO: remove warning
    separate(signature, into = "signature_ct", sep = "[.]", remove = FALSE)

  # Remove signatures which the ctoi score is not the max_score
  sigs_to_remove <- scores_mat_tidy %>%
    group_by(signature_ct, signature) %>%
    mutate(max_score = max(score)) %>%
    filter(signature_ct == mixture_ct & score != max_score)

  scores_mat_tidy <- scores_mat_tidy %>%
    filter(!signature %in% sigs_to_remove$signature) %>%
    select(-mixture_fraction)

  # # Rank by outlier grubbs test
  # scores_mat_tidy %>%
  #   filter(signature_ct == "CL:0000057") %>%
  #   group_by(signature_ct, signature) %>%
  #   summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$p.value) -> y
  #
  # x <- scores_mat_tidy %>%
  #   filter(signature == "CL:0000057.0.585.50%,50%.41") %>%
  #   arrange(-score)
  # x$score
  # res <- outliers::grubbs.test(x$score, type = 10, opposite = FALSE, two.sided = FALSE)
  # res$alternative
  # res$p.value
  # as.numeric(sapply(strsplit(res$alternativ, " "), "[[", 3))


  scores_mat_tidy_ranked <- scores_mat_tidy %>%
    # Rank by CTOI
    filter(signature_ct == mixture_ct) %>%
    group_by(signature_ct) %>%
    mutate(CTOI_rank = percent_rank(score)*100) %>%
    full_join(
      # Rank by CT
      scores_mat_tidy %>%
        filter(signature_ct != mixture_ct) %>%
        group_by(signature_ct, signature) %>%
        summarise(median_score = median(score)) %>%
        mutate(CT_rank = percent_rank(dplyr::desc(median_score))*100) %>%
        select(signature_ct, signature, CT_rank), by = c("signature_ct", "signature")
    ) %>%
    full_join(
      # Rank by diff
      scores_mat_tidy %>%
        mutate(is_ctoi = ifelse(mixture_ct == signature_ct, "yes", "no")) %>%
        group_by(signature_ct, signature, is_ctoi) %>%
        top_n(1, wt=score) %>%
        group_by(signature_ct, signature) %>%
        mutate(top_diff = lag(score) - score) %>%
        drop_na() %>%
        group_by(signature_ct) %>%
        mutate(diff_rank = percent_rank(top_diff)*100) %>%
        select(signature_ct, signature, diff_rank), by = c("signature_ct", "signature")
    ) %>%
    full_join(
      # Rank by outlier grubbs test
      scores_mat_tidy %>%
        group_by(signature_ct, signature) %>%
        summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$p.value) %>%
        group_by(signature_ct) %>%
        mutate(grubbs_rank = percent_rank(dplyr::desc(grubbs_pvalue))*100) %>%
        select(signature_ct, signature, grubbs_rank), by = c("signature_ct", "signature")
    ) %>%
    ungroup()

  return(scores_mat_tidy_ranked)
}




# Testing -----------------------------------



# # x <- scores_mat_tidy %>%
# #   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
# #   group_by(signature_ct, fraction, signature, is_ctoi) %>%
# #   top_n(3, wt=score)
# # x %>%
# #   filter(signature_ct == "CL:0000840", fraction == 1, signature == "CL:0000840.3.25%,75%.60")
#
#
# # Heatmap from score_mat
# scores_mat <- as.matrix(scores_mat)
# scores_mat_sub <- scores_mat[startsWith(rownames(scores_mat), "CL:0000840"), grep("*_1_", colnames(scores_mat))]
# colnames(scores_mat_sub) <- unlist(lapply(strsplit(colnames(scores_mat_sub), "_"), "[", 1))
#
# types2use <- names(sort(celltype_cor_mat["CL:0000840",], decreasing = T))
# types2use <- types2use[!types2use %in% dep_list[["CL:0000840"]]]
# sigsSorted <- scores_mat_tidy_ranked[scores_mat_tidy_ranked$signature_ct== "CL:0000840",]$signature
# scores_mat_sub <- scores_mat_sub[sigsSorted, types2use]
#
# pheatmap::pheatmap(scores_mat_sub, cluster_rows=F, cluster_cols=F)

