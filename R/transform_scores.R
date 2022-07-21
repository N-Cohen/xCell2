


transformScores <- function(scores_mat_frac_tidy){

  # TODO: check mixture_fraction hard numbers in the filters

  # Find calibration values
  calibration_values <- scores_mat_frac_tidy %>%
    filter(signature_ct == mixture_ct) %>%
    # Take the mean of all signatures' scores per cell type per mixture fraction
    group_by(signature_ct, mixture_fraction) %>%
    summarise(mean_score = mean(score)) %>%
    # The shift value is the minimum score per cell type
    mutate(shift_value = min(mean_score)) %>%
    # The scale factor is the slope of the new linear line (i.e., 0.25/shifted_score(25%))
    mutate(shifted_score = mean_score - shift_value) %>%
    mutate(scale_factor = max(as.numeric(mixture_fraction))/max(shifted_score)) %>%
    dplyr::select(signature_ct, shift_value, scale_factor) %>%
    unique()

  # Transform scores
  scores_transformed <- scores_mat_frac_tidy %>%
    filter(mixture_fraction == 0.25) %>%
    group_by(signature_ct, mixture_ct) %>%
    summarise(mean_score = mean(score)) %>%
    left_join(calibration_values, by = "signature_ct") %>%
    rowwise() %>%
    mutate(transformed_score = (mean_score-shift_value)*scale_factor)

  return(list(calibration_values, scores_transformed))

}

#
# # Spillover matrix
# # rows are the scores of cell types based on their signatures
# # cols are the score of cell types based on the signatures of the rows
#
# spillover_mat <- scores_transformed %>%
#   # Transform to wide matrix
#   select(signature_ct, mixture_ct, transformed_score) %>%
#   arrange(signature_ct, mixture_ct) %>%
#   pivot_wider(names_from = mixture_ct, values_from = transformed_score) %>%
#   column_to_rownames(var = "signature_ct") %>%
#   as.matrix()
#
# spillover_mat <- apply(spillover_mat, 2, function(x) {x-median(x, na.rm = TRUE)})
# spillover_mat <- spillover_mat / diag(spillover_mat)
#
# # non-diagonal values greater than 0.5 are corrected to 0.5 (so that spillover correction isn't too strong):
# spillover_mat[spillover_mat>0.5] <- 0.5; diag(spillover_mat) <- 1
# spillover_mat[spillover_mat < 0] <- 0
#
# pheatmap::pheatmap(spillover_mat, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "yellow", "orange", "red"))(75))
#
#
# # Plots
# ctoi <- unique(labels$label)[4]
#
# scores_mat_tidy %>%
#   filter(signature_ct ==  ctoi & signature %in% pull(signatures_filtered, signature) & mixture_fraction < 1) %>%
#   group_by(mixture_ct, mixture_fraction) %>%
#   summarise(mean_score = mean(score)) %>%
#   mutate(R = cor(as.numeric(mixture_fraction), mean_score)) %>%
#   mutate(is_CTOI = ifelse(mixture_ct == ctoi, "yes", "no")) %>%
#   ggplot(., aes(x=mixture_fraction, y=mean_score, col=R, shape=is_CTOI)) +
#   geom_point(size = 3) +
#   scale_color_gradientn(colours = RColorBrewer::brewer.pal(5, "RdBu")) +
#   theme(panel.grid.major = element_line(colour = "gray86"),
#     axis.title = element_text(size = 15),
#     axis.text = element_text(size = 15, colour = "black"),
#     axis.text.x = element_text(size = 12,
#         vjust = 0.5, angle = 40), plot.title = element_text(size = 15),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(size = 15),
#     panel.background = element_rect(fill = NA),
#     legend.key = element_rect(size = 0.9),
#     legend.background = element_rect(fill = NA),
#     legend.position = "bottom", legend.direction = "horizontal") +labs(title = paste0("Mixture scores for ", ctoi), x = "Fraction",
#     y = "Score", shape = "CTOI?")
