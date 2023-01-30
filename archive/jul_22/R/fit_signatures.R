

fitSignatures <- function(signatures_filtered, scores_mat){

  # Make score mat tidy
  scores_mat <- as_tibble(scores_mat, rownames = NA)

  scores_mat_tidy <- scores_mat %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature,
                 names_to = c("mixture_ct", "mixture_fraction", "mixture_control"),
                 names_sep = "_",
                 values_to = "score") %>%
    filter(signature %in% signatures_filtered$signature & mixture_fraction != 1) %>%
    separate(signature, into = "signature_ct", sep = "[.]", remove = FALSE)

  # # Run Lasso regression - canceled (for now)
  # # TODO (!): Check the warnings
  # signatures_fit <- scores_mat_tidy %>%
  #   # filter(signature_ct == mixture_ct) %>%
  #   unite(mixture, mixture_ct:mixture_control, sep = "_") %>%
  #   nest_by(signature_ct) %>%
  #   rename(train = data) %>%
  #   # Make a mixture x signature matrix for each cell type
  #   mutate(train = list(drop_na(pivot_wider(train, names_from = signature, values_from = score)))) %>%
  #   # Use the first column of train matrix (mixture name) as a target vector
  #   mutate(mixtures = list(pull(train[,1])), train = list(train[,-1])) %>%
  #   mutate(train = list(as.matrix(train))) %>%
  #   # Modify target vector: all values are zero expect the CTOI abundance
  #   mutate(traget = list(ifelse(startsWith(mixtures, paste0(signature_ct, "_")), mixtures, "_0_"))) %>%
  #   mutate(traget = list(as.numeric(unlist(lapply(strsplit(traget, "_"), "[", 2))))) %>%
  #   # Run Lasso regression
  #   mutate(lasso_fit = list(glmnet::cv.glmnet(train, traget, type.measure = "mse", alpha = 1))) %>%
  #   # Use Lasso model to run prediction
  #   mutate(predicted_scores = list(as.numeric(predict(lasso_fit, newx = train, s = "lambda.min")))) %>%
  #   # Use lm for scores transformation
  #   mutate(lm = list(lm(traget ~ predicted_scores))) %>%
  #   # ---- TODO
  #   # Make the lowest prediction score represent 0%
  #   # ---
  #   # Transform predicted scores to abundance
  #   mutate(transformed_predictions = list(lm$coefficients[2] * predicted_scores + lm$coefficients[1]))

  # Take the mean of all signatures' scores per cell type per mixture fraction
  mean_scores_mat_tidy <- scores_mat_tidy %>%
    filter(signature_ct == mixture_ct) %>%
    group_by(signature_ct, mixture_fraction) %>%
    summarise(mean_score = mean(score)) %>%
    ungroup()

  # Get calibration parameters
  calibration_data <- mean_scores_mat_tidy %>%
    # The shift value is the minimum score per cell type
    group_by(signature_ct) %>%
    mutate(shift_value = min(mean_score, na.rm = TRUE)) %>%
    select(signature_ct, shift_value) %>%
    unique() %>%
    ungroup()

  calibration_data <- mean_scores_mat_tidy %>%
      filter(mixture_fraction == max(scores_mat_tidy$mixture_fraction)) %>%
      # The scale factor is the slope of the new linear line (i.e., 0.25/shifted_score(25%))
      left_join(calibration_data, by = "signature_ct") %>%
      mutate(shifted_score = mean_score - shift_value) %>%
      mutate(scale_factor = as.numeric(mixture_fraction)/shifted_score) %>%
      select(signature_ct, shift_value, scale_factor)


 return(calibration_data)
}


# Plot before/after calibration
mean_scores_mat_tidy %>%
  filter(mixture_fraction >= 0.01) %>%
  left_join(calibration_data, by="signature_ct") %>%
  mutate(shifted_score = mean_score - shift_value) %>%
  mutate(adjusted_score = scale_factor*shifted_score) %>%
  select(signature_ct, mixture_fraction, mean_score, shifted_score, adjusted_score) %>%
  filter(signature_ct == "CL:0000066") %>%
  pivot_longer(cols = mean_score:adjusted_score, names_to = "transformation", values_to = "score") %>%
  ggplot(., aes(x=as.numeric(mixture_fraction), y=score, col=transformation)) +
  geom_point()


# Plot linear dependencies between cell types
scores_mat_tidy %>%
  filter(signature_ct == "CL:0000863" & mixture_ct %in% c("CL:0000863", "CL:0000890")) %>%
  group_by(signature_ct, mixture_ct, mixture_fraction) %>%
  summarise(mean_score = mean(score)) %>%
  ungroup() %>%
  filter(mixture_fraction >= 0.01) %>%
  left_join(calibration_data, by="signature_ct") %>%
  mutate(adjusted_score = (mean_score-shift_value)*scale_factor) %>%
  select(mixture_ct, mixture_fraction, adjusted_score) %>%
  pivot_wider(names_from = mixture_ct, values_from = adjusted_score) -> x

colnames(x)[2:3] <- c("M1", "M2")

  ggplot(x, aes(x=M1, y=M2)) +
  geom_point()






# # Plot linear models
# scores_mat_tidy %>%
#   filter(signature_ct == mixture_ct) %>%
#   group_by(signature_ct, mixture_fraction) %>%
#   mutate(mean_score = mean(score)) %>%
#   select(-signature, -mixture_ct, -mixture_control, -score) %>%
#   group_by(signature_ct) %>%
#   mutate(shift_value = min(mean_score)) %>%
#   mutate(score_shifted = mean_score - shift_value) %>%
#   select(-mean_score) %>%
#   group_by(signature_ct, shift_value) %>%
#   # Choose cell type
#   filter(signature_ct == "CL:0000624") %>%
#   ggplot(., aes(x = score_shifted, y = mixture_fraction)) +
#   geom_point() +
#   stat_smooth(method = "lm", col = "red")
#



# # Plot score to abundance
# ct <- "CL:0000624"
#
# scores_mat_tidy %>%
#   filter(signature_ct == ct) %>%
#   mutate(is_CTOI = ifelse(signature_ct == mixture_ct, "yes", "no")) %>%
#   mutate(mixture_fraction= ifelse(signature_ct == mixture_ct, mixture_fraction, 0)) %>%
#   select(score, mixture_fraction, is_CTOI) %>%
#   ggplot(., aes(x=score, y=mixture_fraction, col=is_CTOI)) +
#   geom_point()
#
# scores_mat_tidy %>%
#   filter(signature_ct == ct & signature_ct == mixture_ct & signature == "CL:0000815.1.10%,90%.41") %>%
#   select(score, mixture_fraction) %>%
#   ggplot(., aes(x=score, y=mixture_fraction)) +
#   geom_point()
#
#
#
# # # Check linear transformation
# ps <- signatures_fit[1,]$predicted_scores[[1]]
# tp <- signatures_fit[1,]$transformed_predictions[[1]]
# is_ctoi <- ifelse(signatures_fit[1,]$traget[[1]] > 0, "yes", "no")
#
# df <- data.frame(ps = ps, tp = tp, is_ctoi = is_ctoi)
#
# library(ggThemeAssist)
# ggplot(df, aes(ps, tp, col=is_ctoi)) +
#   geom_point() +
#   scale_x_continuous(breaks=seq(-25, 0.25, 1))
#   theme(axis.text = element_text(size = 13),
#     panel.background = element_rect(fill = "gray97"),
#     plot.background = element_rect(fill = "white")) +labs(x = "Predicted Scores", y = "Transformed Scores",
#     colour = "CTOI")
#



# # With tidymodels
# library(tidymodels)
#
# lasso_model <- linear_reg(penalty = 0.01, mixture = 1) %>% # lasso: 1, ridge: 0
#   set_engine("glmnet")
#
# x <- scores_mat_tidy %>%
#   mutate(ctoi_fraction =  as.numeric(ifelse(signature_ct == mixture_ct, mixture_fraction, 0))) %>%
#   select(signature_ct, signature, mixture_ct, ctoi_fraction, score) %>%
#   pivot_longer(cols = ctoi_fraction:score) %>%
#   mutate(signature = ifelse(signature_ct == mixture_ct, name, signature)) %>%
#   select(-mixture_ct, -name) %>%
#   nest_by(signature_ct) %>%
#   mutate(data = list(drop_na(pivot_wider(data, names_from = signature, values_from = score))))
#



# Heatmap
# chnageLabel <- function(ref, ont, label = "fine"){
#   if (label == "main") {
#     out <- unique(ref$label.main[ref$label.ont == ont])
#     return(out[!is.na(out)][1])
#   }else if(label == "fine"){
#     out <- unique(ref$label.fine[ref$label.ont == ont])
#     return(out[!is.na(out)][1])
#   }
# }
#
#
# plots <- list()
#
# x1 <- signatures_fit %>%
#   select(signature_ct, mixtures, predicted_scores) %>%
#   unnest(cols = c(mixtures, predicted_scores)) %>%
#   pivot_wider(names_from = signature_ct, values_from = predicted_scores) %>%
#   select(-mixtures) %>%
#   select(celltypes)
#
# colnames(x1) <- unlist(lapply(colnames(x1), function(x) chnageLabel(ref, x)))
# plots[["predicted"]] <- pheatmap::pheatmap(x1, cluster_rows=F, cluster_cols=F, show_rownames=F)[[4]]
#
# x2 <- signatures_fit %>%
#   select(signature_ct, mixtures, transformed_predictions) %>%
#   unnest(cols = c(mixtures, transformed_predictions)) %>%
#   pivot_wider(names_from = signature_ct, values_from = transformed_predictions) %>%
#   select(-mixtures) %>%
#   select(celltypes)
#
# colnames(x2) <- unlist(lapply(colnames(x2), function(x) chnageLabel(ref, x)))
# plots[["transformed"]] <- pheatmap::pheatmap(x2, cluster_rows=F, cluster_cols=F, show_rownames=F)[[4]]
#
# gridExtra::grid.arrange(grobs=plots, ncol=2)
