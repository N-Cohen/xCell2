spillOverMat <- function(signatures_fit, mixture_fractions, celltype_cor_mat = celltype_cor_mat, cor_cells_cutoff = .92){


  # Make score mat tidy
  scores_mat <- as_tibble(scores_mat, rownames = NA)

  scores_mat_tidy <- scores_mat %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature,
                 names_to = c("mixture_ct", "mixture_fraction", "mixture_control"),
                 names_sep = "_",
                 values_to = "score") %>%
    filter(signature %in% signatures_filtered$signature & mixture_fraction == 0.25) %>%
    separate(signature, into = "signature_ct", sep = "[.]", remove = FALSE)


  scores_mat_tidy %>%
    unite(mixture, mixture_ct:mixture_control, sep = "_") %>%
    # Nest scores by signature ct
    nest_by(signature_ct) %>%
    rename(new_x = data) %>%
    # new_x is a wide format matrix of signatures (columns) x 25% CTOI mixtures (rows)
    mutate(new_x = list(drop_na(pivot_wider(new_x, names_from = signature, values_from = score)))) %>%
    # Pull the first column (mixtures ID) from the score matrix
    mutate(mixtures = list(pull(new_x[,1])), new_x = list(new_x[,-1])) %>%
    mutate(new_x = list(as.matrix(new_x))) %>%
    # Join lasso models for each CT
    left_join(signatures_fit %>% select(lasso_fit), by = "signature_ct") %>%
    # Run predictions with the lasso model
    mutate(all_ct_prediction = list(as.numeric(predict(lasso_fit, newx = new_x, s = "lambda.min")))) %>%
    select(signature_ct, mixtures,all_ct_prediction) %>%
    unnest(cols = c(mixtures, all_ct_prediction)) %>%
    separate(mixtures, into = "mixture_ct", sep = "_", remove = FALSE) %>%
    select(signature_ct, mixture_ct, all_ct_prediction) %>%
    arrange(signature_ct, mixture_ct) %>%
    pivot_wider(names_from = mixture_ct, values_from = all_ct_prediction) -> x





  highest_frac_index <- which.max(mixture_fractions)

  split_mat <- signatures_fit %>%
    mutate(transformed_predictions = list(transformed_predictions[seq(highest_frac_index, length(transformed_predictions), highest_frac_index)]),
           mixture_ct = list(unique(unlist(lapply(mixtures, function(x) sub("_.*", "", x)))))) %>%
    select(signature_ct, mixture_ct, transformed_predictions) %>%
    unnest(cols = c(mixture_ct, transformed_predictions)) %>%
    arrange(signature_ct, mixture_ct) %>%
    pivot_wider(names_from = mixture_ct, values_from = transformed_predictions)

  split_mat.df <- as.data.frame(split_mat)
  rownames(split_mat.df) <- split_mat.df$signature_ct
  split_mat.df <- split_mat.df[,-1]

  split_mat <- as.matrix(split_mat.df)
  split_mat <- split_mat / diag(split_mat)


  split_mat[split_mat > 0.5] <- 0.5
  split_mat[split_mat < 0] <- 0
  diag(split_mat) <- 1

  cor_mat <- celltype_cor_mat[rownames(split_mat), colnames(split_mat)] < cor_cells_cutoff
  cor_mat[upper.tri(cor_mat, diag = FALSE)] <- FALSE
  split_mat[cor_mat] <- 0


}

chnageLabel <- function(ref, ont, label = "fine"){
  if (label == "main") {
    out <- unique(ref$label.main[ref$label.ont == ont])
    return(out[!is.na(out)][1])
  }else if(label == "fine"){
    out <- unique(ref$label.fine[ref$label.ont == ont])
    return(out[!is.na(out)][1])
  }
}

colnames(split_mat) <- unlist(lapply(colnames(split_mat), function(x) chnageLabel(ref, x)))
rownames(split_mat) <- colnames(split_mat)

pheatmap::pheatmap(split_mat, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "yellow", "orange", "red"))(75))

