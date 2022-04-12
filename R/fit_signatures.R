

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

  # Run Lasso regression
  # TODO (!): Check the warnings
  signatures_fit <- scores_mat_tidy %>%
    filter(signature_ct == mixture_ct) %>%
    unite(mixture, mixture_ct:mixture_control, sep = "_") %>%
    nest_by(signature_ct) %>%
    rename(train = data) %>%
    # Make a mixture x signature matrix for each cell type
    mutate(train = list(drop_na(pivot_wider(train, names_from = signature, values_from = score)))) %>%
    # Use the first column of train matrix (mixture name) as a target vector
    mutate(mixtures = list(pull(train[,1])), train = list(train[,-1])) %>%
    mutate(train = list(as.matrix(train))) %>%
    # Modify target vector: all values are zero expect the CTOI abundance
    # Notice (!): We are only using mixtures of CTOI, therefore no need for the "zeros"
    mutate(traget = list(ifelse(startsWith(mixtures, paste0(signature_ct, "_")), mixtures, "_0_"))) %>%
    mutate(traget = list(as.numeric(unlist(lapply(strsplit(traget, "_"), "[", 2))))) %>%
    # Run Lasso regression
    mutate(lasso_fit = list(glmnet::cv.glmnet(train, traget, type.measure = "mse", alpha = 1))) %>%
    # Use Lasso model to run prediction
    mutate(predicted_scores = list(as.numeric(predict(lasso_fit, newx = train, s = "lambda.min")))) %>%
    # Use lm for scores transformation
    mutate(lm = list(lm(traget ~ predicted_scores))) %>%
    # Transform predicted scores to abundance
    mutate(transformed_predictions = list(lm$coefficients[2] * predicted_scores + lm$coefficients[1]))

 return(signatures_fit)
}


# # Check linear transformation
ps <- signatures_fit[1,]$predicted_scores[[1]]
tp <- signatures_fit[1,]$transformed_predictions[[1]]
is_ctoi <- ifelse(signatures_fit[1,]$traget[[1]] > 0, "yes", "no")

df <- data.frame(ps = ps, tp = tp, is_ctoi = is_ctoi)

library(ggThemeAssist)
ggplot(df, aes(ps, tp, col=is_ctoi)) +
  geom_point() +
  scale_x_continuous(breaks=seq(-25, 0.25, 1))
  theme(axis.text = element_text(size = 13),
    panel.background = element_rect(fill = "gray97"),
    plot.background = element_rect(fill = "white")) +labs(x = "Predicted Scores", y = "Transformed Scores",
    colour = "CTOI")




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
