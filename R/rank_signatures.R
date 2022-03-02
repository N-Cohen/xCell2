library(tidyverse)


rankSignatures <- function(scores_mat){

  scores_mat <- as_tibble(scores_mat, rownames = NA)

  scores_mat_tidy <- scores_mat %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature,
                 names_to = c("mixture_ctoi", "fraction", "control"),
                 names_sep = "_",
                 values_to = "score") %>%
    filter(fraction == 1)

  # TODO: make this run faster
  scores_mat_tidy$signature_ct <- unlist(lapply(strsplit(scores_mat_tidy$signature, "[.]"), "[", 1))


  scores_mat_tidy_ranked <- scores_mat_tidy %>%
    # Rank by CTOI
    filter(signature_ct == mixture_ctoi) %>%
    group_by(signature_ct) %>%
    mutate(CTOI_rank = percent_rank(score)*100) %>%
    arrange(-CTOI_rank) %>%
    full_join(
      # Rank by CT
      scores_mat_tidy %>%
        filter(signature_ct != mixture_ctoi) %>%
        drop_na() %>%
        group_by(signature_ct, signature) %>%
        summarise(median_score = median(score)) %>%
        mutate(CT_rank = percent_rank(desc(median_score))*100) %>%
        arrange(-CT_rank), by = c("signature_ct", "signature")
    ) %>%
    full_join(
      # Rank by diff
      scores_mat_tidy %>%
        filter(fraction == 1) %>%
        mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
        group_by(signature_ct, signature, is_ctoi) %>%
        top_n(1, wt=score) %>%
        group_by(signature_ct, signature) %>%
        mutate(top_diff = lag(score) - score) %>%
        drop_na() %>%
        group_by(signature_ct) %>%
        mutate(diff_rank = percent_rank(top_diff)*100) %>%
        arrange(-diff_rank), by = c("signature_ct", "signature")
    ) %>%
    ungroup()

  return(scores_mat_tidy_ranked)
}



# # Rank by CTOI
# scores_mat_tidy %>%
#   filter(signature_ct == mixture_ctoi & fraction == 1) %>%
#   group_by(signature_ct) %>%
#   mutate(CTOI_rank = percent_rank(score)*100) %>%
#   arrange(-CTOI_rank)
#
# # Rank by CT
# scores_mat_tidy %>%
#   filter(signature_ct != mixture_ctoi & fraction == 1) %>%
#   drop_na() %>%
#   group_by(signature_ct, signature) %>%
#   summarise(median_score = median(score)) %>%
#   mutate(CT_rank = percent_rank(desc(median_score))*100) %>%
#   arrange(-CT_rank)
#
# # Rank by diff
# scores_mat_tidy %>%
#   filter(fraction == 1) %>%
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
#   group_by(signature_ct, signature, is_ctoi) %>%
#   top_n(1, wt=score) %>%
#   group_by(signature_ct, signature) %>%
#   mutate(top_diff = lag(score) - score) %>%
#   drop_na() %>%
#   group_by(signature_ct) %>%
#   mutate(diff_rank = percent_rank(top_diff)*100) %>%
#   arrange(-diff_rank)


# Testing -----------------------------------
# scores_mat <- as_tibble(scores_mat, rownames = NA)
#
# scores_mat_tidy <- scores_mat %>%
#   rownames_to_column(var = "signature") %>%
#   pivot_longer(cols = -signature,
#                names_to = c("mixture_ctoi", "fraction", "control"),
#                names_sep = "_",
#                values_to = "score")
#
# scores_mat_tidy$signature_ct <- unlist(lapply(strsplit(scores_mat_tidy$signature, "[.]"), "[", 1))
#
#
# # Ranking:
# # (1) CTOI_rank: for each cell type, fraction rank signatures by highest CTOI score then sum
# scores_mat_tidy %>%
#   filter(signature_ct == mixture_ctoi) %>%
#   group_by(signature_ct, fraction) %>%
#   mutate(rank = dense_rank(score)) %>%
#   group_by(signature_ct, signature) %>%
#   summarise(CTOI_rank = sum(rank)) %>%
#   arrange(-CTOI_rank)
#
# # (2) CT_rank: for each cell type, fraction, rank signatures by lowest median CT score then sum
# # TODO: CL:0000017 is missing!
# scores_mat_tidy %>%
#   filter(signature_ct != mixture_ctoi) %>%
#   group_by(signature_ct, fraction, signature) %>%
#   summarise(median_score = median(score)) %>%
#   drop_na() %>%
#   mutate(rank = dense_rank(desc(median_score))) %>%
#   group_by(signature_ct, signature) %>%
#   summarise(CT_rank = sum(rank)) %>%
#   arrange(-CT_rank)
#
#
# # (3) diff_rank: for each cell type, fraction, rank signatures by highest difference between scores of CTOI to the next highest CT
# scores_mat_tidy %>%
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
#   group_by(signature_ct, fraction, signature, is_ctoi) %>%
#   top_n(1, wt=score) %>%
#   group_by(signature_ct, fraction, signature) %>%
#   mutate(top_diff = lag(score) - score) %>%
#   drop_na() %>%
#   group_by(signature_ct, fraction) %>%
#   mutate(rank = dense_rank(top_diff)) %>%
#   group_by(signature_ct, signature) %>%
#   summarise(diff_rank = sum(rank)) %>%
#   arrange(-diff_rank)
#
#
#
#
#
#
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

