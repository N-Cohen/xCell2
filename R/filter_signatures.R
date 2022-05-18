

filterSignatures <- function(scores_mat_tidy, take_top_per){

  # Rank signature by Grubbs' test
  signatures_filtered <- scores_mat_tidy %>%
    filter(mixture_fraction == 1) %>%
    drop_na() %>%
    # Remove signatures which the ctoi score is not the max_score
    group_by(signature_ct, signature) %>%
    filter(any(signature_ct == mixture_ct & score == max(score))) %>%
    # Rank signatures
    summarise(grubbs_pvalue = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$p.value) %>%
    mutate(grubbs_rank = percent_rank(dplyr::desc(grubbs_pvalue))*100) %>%
    arrange(-grubbs_rank, .by_group = TRUE) %>%
    # Filter signatures
    filter(grubbs_rank >= quantile(grubbs_rank, 1-take_top_per))

  return(signatures_filtered)
}




# old
# ranks_weights = c(0, 0, 0, 1)
# signatures_filtered <- signatures_ranked %>%
#   mutate(CT_rank = coalesce(CT_rank, CTOI_rank), diff_rank = coalesce(diff_rank, CTOI_rank)) %>%
#   mutate(final_rank = ranks_weights[1]*CTOI_rank + ranks_weights[2]*CT_rank + ranks_weights[3]*diff_rank + ranks_weights[4]*grubbs_rank) %>%
#   group_by(signature_ct) %>%
#   arrange(-final_rank, .by_group = TRUE) %>%
#   filter(final_rank > quantile(final_rank, 1-take_top_per)) %>%
#   select(signature_ct, signature) %>%
#   ungroup()
















#
#
# library(tidyverse)
#
#
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
# # Check correlation between fraction and score
# scores_mat_tidy %>%
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, TRUE, FALSE)) %>%
#   filter(signature == "CL:0002328.5.25%,75%.62") %>%
#   mutate(fraction = ifelse(mixture_ctoi == signature_ct, fraction, 0)) %>%
#   ggplot(., aes(x=score, y=fraction, col=is_ctoi)) +
#   geom_point()
#
# # Five signatures have negative correlation (?)
# scores_mat_tidy %>%
#   filter(mixture_ctoi == signature_ct & fraction != 1) %>%
#   group_by(signature) %>%
#   mutate(fraction = as.numeric(fraction)) %>%
#   summarise(cor = cor(fraction, score)) %>%
#   arrange(cor)
#
# scores_mat_tidy %>%
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, TRUE, FALSE)) %>%
#   filter(signature == "CL:0000576.0.25%,75%.61") %>%
#   mutate(fraction = ifelse(mixture_ctoi == signature_ct, fraction, 0)) %>%
#   ggplot(., aes(x=score, y=fraction, col=is_ctoi)) +
#   geom_point()
#
# signatures_list$`CL:0000625`$`0`$`33%,67%`$`59`
# x <- scores_mat_tidy %>%
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, TRUE, FALSE)) %>%
#   filter(signature == "CL:0000625.0.33%,67%.59")
#
#
# # Rule #1 - CTOI score must be the highest:
# x <- scores_mat_tidy %>%
#   group_by(signature, fraction) %>%
#   top_n(1, wt=score) %>%
#   mutate(ctoi_is_top = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
#   filter(ctoi_is_top == "no")
#
# unique(x$signature) # Remove those signatures
# scores_mat_tidy_rule_one <- scores_mat_tidy %>%
#   filter(!signature %in% unique(x$signature))
#
# # Check how many signatures left
# x <- scores_mat_tidy_rule_one %>%
#   select(signature_ct, signature) %>%
#   unique() %>%
#   group_by(signature_ct) %>%
#   summarise(n = n())
#
# celltypes[!celltypes %in% x$signature_ct] # 8 cell types lost all signatures
#
# # Rule #2 - Score of lowest fraction of CTOI must be higher than max(score) of any other CT
# x <- scores_mat_tidy %>%
#   filter(signature_ct != "CL:0000017") %>% # This CT is a problem because he is control for all others
#   mutate(is_ctoi = ifelse(mixture_ctoi == signature_ct, "yes", "no")) %>%
#   group_by(signature, is_ctoi) %>%
#   slice(case_when(any(is_ctoi == "yes") ~ which.min(score),
#                   any(is_ctoi == "no") ~ which.max(score),
#                   TRUE ~ which.min(score))) %>%
#   group_by(signature) %>%
#   mutate(diff = score - lag(score)) %>%
#   drop_na()
#
#
#
#
#
# # Rule #2 - Top diff score
# scores_mat_tidy_rule_one %>%
#   filter(fraction == 1) %>%
#   group_by(signature) %>%
#   top_n(2, wt=score) %>%
#   arrange(-score, by_group = TRUE) %>%
#   mutate(top_diff = lag(score) - score) %>%
#   drop_na() %>%
#   ungroup() %>%
#   ggplot(., aes(x=signature_ct, y=top_diff)) +
#   geom_boxplot()
#
#
#
#
#
#
#
#
#
#
#
# # Old ------------------------------------------
#
# library(tidyverse)
# library(singscore)
# library(pheatmap)
#
# signatures_list <- readRDS("data_for_dev/signatures_list.RData")
# ref <- readRDS("~/Documents/xCell2.0/HPCA.RData")
# OBOfile <- "~/Documents/xCell2.0/cl.obo"
# counts <- ref@assays@data$logcount
# samples <- ref$label.ont
# celltypes <- unique(samples[!is.na(samples)])
#
# dep_list <- readRDS("data_for_dev/dep_list.RData")
# celltype_cor_mat <- readRDS("data_for_dev/celltype_cor_mat.RData")
#
#
# # This function return a table of signatures and their scores
# scoreSignatures <- function(counts, samples, celltypes, signatures_list, dep_list, celltype_cor_mat){
#
#   ctoi_vec <- c()
#   target_celltype_vec <- c()
#   signature_vec <- c()
#   median_scores_vec <- c()
#
#   for (ctoi in names(signatures_list)) {
#     message(paste0("Calculating scores for ", ctoi))
#     dep_cells <- dep_list[[ctoi]]
#     celltypes_nodep <- celltypes[!celltypes %in% dep_cells]
#     all_sigs <- lapply(rapply(signatures_list[[ctoi]], enquote, how="unlist"), eval)
#
#
#     for (i in 1:length(celltypes_nodep)) {
#
#       # Make separate ranked counts matrix for each cell type
#       ct <- celltypes_nodep[i]
#       counts_ct <- counts[,samples == ct]
#       counts_ct <- counts_ct[,!is.na(colnames(counts_ct))]
#
#       if (sum(samples == ct, na.rm = TRUE) == 1) {
#         counts_ct <- cbind(counts_ct, counts_ct)
#         colnames(counts_ct) <- c("a", "b")
#       }
#       counts_ct_ranked <- rankGenes(counts_ct) # Rank counts
#
#       # Get the median score for each signature
#       median_scores_tmp <- c()
#       for (j in 1:length(all_sigs)) {
#         scores_out <- simpleScore(counts_ct_ranked, upSet = all_sigs[j][[1]])$TotalScore
#         median_scores_tmp <- c(median_scores_tmp, median(scores_out))
#       }
#
#       # Save results
#       ctoi_vec <- c(ctoi_vec, rep(ctoi, length(all_sigs)))
#       target_celltype_vec <- c(target_celltype_vec, rep(ct, length(all_sigs)))
#       signature_vec <- c(signature_vec, names(all_sigs))
#       median_scores_vec <- c(median_scores_vec, median_scores_tmp)
#     }
#
#   }
#
#   score_out <- tibble(CTOI = ctoi_vec, celltype = target_celltype_vec, signature = signature_vec,
#                          median_score = median_scores_vec)
#
# }
#
# # This function return a list of filtered signatures for each cell type
# filterSignatures <- function(sig_scores, signatures_list, take_top = .9){
#   celltypes <- unique(sig_scores$CTOI)
#   filtered_signatures_out <- vector(mode = "list", length = length(celltypes))
#   names(filtered_signatures_out) <- celltypes
#
#   for (ctoi in celltypes) {
#
#     ctoi_df <- sig_scores[sig_scores$CTOI == ctoi,]
#
#     filtered_ctoi_df <- ctoi_df %>%
#       # Rank by CTOI
#       filter(CTOI == celltype) %>%
#       mutate(CTOI_rank =  dense_rank(median_score)) %>%
#       full_join(
#         # Rank by CT
#         ctoi_df %>%
#           filter(CTOI != celltype) %>%
#           group_by(signature) %>%
#           summarise(ct_median_score = median(median_score)) %>%
#           mutate(CT_rank =  dense_rank(desc(ct_median_score))), by = "signature"
#       ) %>%
#       full_join(
#         # Rank by diff
#         ctoi_df %>%
#           group_by(signature) %>%
#           top_n(2, wt=median_score) %>%
#           mutate(top_diff = lag(median_score) - median_score) %>%
#           drop_na() %>%
#           select(signature, top_diff) %>%
#           ungroup() %>%
#           mutate(diff_rank =  dense_rank(top_diff)), by = "signature"
#       ) %>%
#       # Merge all ranks
#       mutate(final_rank = .5*CTOI_rank + .25*CT_rank + .25*diff_rank) %>%
#       arrange(-final_rank) %>%
#       filter(final_rank > quantile(final_rank, take_top))
#
#     all_sigs <- lapply(rapply(signatures_list[[ctoi]], enquote, how="unlist"), eval)
#     filtered_signatures_out[[ctoi]] <- all_sigs[names(all_sigs) %in% filtered_ctoi_df$signature]
#
#   }
#
#   return(filtered_signatures_out)
# }
#
#
# # Testing --------------------------------------
#
# # sig_scores.df <- scoreSignatures(counts, samples, celltypes, signatures_list, dep_list, celltype_cor_mat)
# sig_scores.df <- readRDS("data_for_dev/sig_scores.df.RData")
#
# ctoi <- celltypes[1]
# sig_scores.df <- sig_scores.df %>%
#   filter(CTOI == ctoi)
#
#
# tmp.df <- sig_scores.df %>%
#   filter(CTOI == celltype) %>%
#   mutate(CTOI_rank =  dense_rank(median_score)) %>%
#   full_join(
#     sig_scores.df %>%
#       filter(CTOI != celltype) %>%
#       group_by(signature) %>%
#       summarise(ct_median_score = median(median_score)) %>%
#       mutate(CT_rank =  dense_rank(desc(ct_median_score))), by = "signature"
#   ) %>%
#   full_join(
#     sig_scores.df %>%
#       group_by(signature) %>%
#       top_n(2, wt=median_score) %>%
#       mutate(top_diff = lag(median_score) - median_score) %>%
#       drop_na() %>%
#       select(signature, top_diff) %>%
#       ungroup() %>%
#       mutate(diff_rank =  dense_rank(top_diff)), by = "signature"
#   ) %>%
#   mutate(final_rank = .5*CTOI_rank + .25*CT_rank + .25*diff_rank) %>%
#   arrange(-final_rank) %>%
#   filter(final_rank > quantile(final_rank, .9))
#
#
# sig_scores.df$signature <- factor(sig_scores.df$signature, levels = tmp.df$signature)
# sig_scores.df$celltype <- factor(sig_scores.df$celltype, levels = names(sort(celltype_cor_mat[ctoi,], decreasing = TRUE)))
#
#
# sig_scores.df %>%
#   filter(signature %in% tmp.df$signature) %>%
#   select(celltype, signature, median_score) %>%
#   spread(key = celltype, value = median_score) %>%
#   column_to_rownames(var = "signature") %>%
#   pheatmap(cluster_rows=F, cluster_cols=F)
#
#
#
#
#
# filtered_signatures_list <- filterSignatures(sig_scores, signatures_list, take_top = .9)
#
#
# # tmp.df <- score_out %>%
# #   group_by(signature) %>%
# #   top_n(2, wt=median_score) %>%
# #   mutate(top_diff = lag(median_score) - median_score) %>%
# #   drop_na() %>%
# #   select(signature, top_diff) %>%
# #   ungroup() %>%
# #   mutate(diff_rank =  dense_rank(top_diff))
# #
# #
# # score_out$signature <- factor(score_out$signature, levels = tmp.df$signature)
# #
# # # (1) Score of CTOI must be in top 75% quantile
# # tmp.df <-score_out %>%
# #   filter(is_ctoi == "yes" & median_score >= quantile(score_out[score_out$is_ctoi == "yes",]$median_score, .75))
# # score_out <- score_out[score_out$signature %in% tmp.df$signature,]
# #
# #
# # # (2) Score of all CTs must be in bottom 75% quantile
# # tmp.df <- score_out %>%
# #   filter(is_ctoi == "no") %>%
# #   group_by(signature) %>%
# #   summarise(all_ct_median = median(median_score))
# # score_out <- score_out[score_out$signature %in% tmp.df[tmp.df$all_ct_median < quantile(tmp.df$all_ct_median, .75),]$signature,]
# #
# #
#
#
#
#
#
#
# # Load data
# HPCAsig <- readRDS("~/Documents/xCell2.0/HPCA_signatures.RData")
# hpca <- readRDS("~/Documents/xCell2.0/HPCA.RData")
# counts <- hpca@assays@data$logcount
# celltype_cor_mat <- readRDS("~/Documents/xCell2.0/celltype_cor_mat.RData")
# celltypes <- unique(hpca$label.ont)[!is.na(unique(hpca$label.ont))]
#
# # For example "CL:0000840":
# ctoi <- "CL:0000840"
# dep_cells <- getDependencies(OBOfile = "~/Documents/xCell2.0/cl.obo", celltypes, type = ctoi)
# celltypes_nodep <- celltypes[!celltypes %in% dep_cells]
# all_sigs <- lapply(rapply(HPCAsig$HPCA[[ctoi]], enquote, how="unlist"), eval)
#
#
# # Get scores for each signatures of "CL:0000840":
# celltypes_vec <- c()
# signature_vec <- c()
# median_scores_vec <- c()
#
# for (i in 1:length(celltypes_nodep)) {
#   # Make seperate ranked counts matrix for each cell type
#   ct <- celltypes_nodep[i]
#   counts_ct <- counts[,hpca$label.ont == ct]
#   counts_ct <- counts_ct[,!is.na(colnames(counts_ct))]
#
#   if (sum(hpca$label.ont == ct, na.rm = TRUE) == 1) {
#     counts_ct <- cbind(counts_ct, counts_ct)
#     colnames(counts_ct) <- c("a", "b")
#   }
#   counts_ct_ranked <- rankGenes(counts_ct) # Rank counts
#
#   # Get the median score and SD for each signature
#   median_scores_tmp <- c()
#   sd_scores_tmp <- c()
#   for (j in 1:length(all_sigs)) {
#     scores_out <- simpleScore(counts_ct_ranked, upSet = all_sigs[j][[1]])$TotalScore
#     median_scores_tmp <- c(median_scores_tmp, median(scores_out))
#     sd_scores_tmp <- c(sd_scores_tmp, sd(scores_out))
#   }
#
#   # Save results
#   celltypes_vec <- c(celltypes_vec, rep(ct, length(all_sigs)))
#   signature_vec <- c(signature_vec, names(all_sigs))
#   median_scores_vec <- c(median_scores_vec, median_scores_tmp)
# }
#
#
# df <- tibble(celltype = celltypes_vec, signature = signature_vec,
#              median_score = median_scores_vec)
#
# scores.df <- as.tibble(as.list(celltype_cor_mat[,ctoi])) %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   as.tibble() %>%
#   rename("celltype"=1, "correlation"=2) %>%
#   right_join(df, by = "celltype") %>%
#   mutate(is_ctoi = ifelse(celltype == ctoi, "yes", "no"))
#
#
# # Rank signatures by:
#
# # (1) Score of CTOI must be in top 75% quartile
# tmp.df <-scores.df %>%
#   filter(is_ctoi == "yes" & median_score >= quantile(scores.df[scores.df$is_ctoi == "yes",]$median_score, .75))
# scores.df <- scores.df[scores.df$signature %in% tmp.df$signature,]
#
#
# # (2) Score of all CTs must be in bottom 75% quartile
# tmp.df <- scores.df %>%
#   filter(is_ctoi == "no") %>%
#   group_by(signature) %>%
#   summarise(all_ct_median = median(median_score))
#
# scores.df <- scores.df[scores.df$signature %in% tmp.df[tmp.df$all_ct_median < quantile(tmp.df$all_ct_median, .75),]$signature,]
#
#
# # (3) Difference between CTOI to the top CT
# tmp.df <- scores.df %>%
#   group_by(signature) %>%
#   top_n(2, wt=median_score) %>%
#   arrange(signature) %>%
#   mutate(Diff = lag(median_score) - median_score) %>%
#   drop_na() %>%
#   arrange(-Diff)
#
#
# scores.df$celltype <- factor(scores.df$celltype, levels = names(sort(celltype_cor_mat[ctoi,], decreasing = TRUE)))
# scores.df$signature <- factor(scores.df$signature, levels = tmp.df$signature)
#
# ggplot(scores.df, aes(x=signature, y=median_score)) +
#   geom_boxplot(aes(col=is_ctoi)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
#
# scores.df %>%
#   select(celltype, signature, median_score) %>%
#   spread(key = celltype, value = median_score) %>%
#   column_to_rownames(var = "signature") %>%
#   pheatmap(cluster_rows=F, cluster_cols=F)
#
#
#
#
#
