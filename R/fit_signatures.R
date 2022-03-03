library(glmnet)


fitSignatures <- function(signatures_filtered, scores_mat){

  # Make score mat tidy
  scores_mat <- as_tibble(scores_mat, rownames = NA)
  scores_mat_tidy <- scores_mat %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature,
                 names_to = c("mixture_ctoi", "fraction", "control"),
                 names_sep = "_",
                 values_to = "score") %>%
    filter(signature %in% signatures_filtered$signature & fraction != 1)
  scores_mat_tidy$signature_ct <- unlist(lapply(strsplit(scores_mat_tidy$signature, "[.]"), "[", 1))

  # Run Lasso regression
  # TODO (!): Check the warnings
  signatures_fitted <- scores_mat_tidy %>%
    unite(mixture, mixture_ctoi:control, sep = "_") %>%
    nest_by(signature_ct) %>%
    rename(train = data) %>%
    # Make a mixture x signature matrix for each cell type
    mutate(train = list(drop_na(pivot_wider(train, names_from = signature, values_from = score)))) %>%
    # Use the first column of train matrix (mixture name) as a target vector
    mutate(traget = list(pull(train[,1])), train = list(train[,-1])) %>%
    mutate(train = list(as.matrix(train))) %>%
    # Modify target vector: all values are zero expect the CTOI abundace
    mutate(traget = list(ifelse(startsWith(traget, paste0(signature_ct, "_")), traget, "_0_"))) %>%
    mutate(traget = list(as.numeric(unlist(lapply(strsplit(traget, "_"), "[", 2))))) %>%
    # Run Lasso regression
    mutate(lasso_fit = list(cv.glmnet(train, traget, type.measure = "mse", alpha = 1)))

 return(signatures_fitted)
}


# set.seed(123)
#
# sigs <- signatures_filtered[signatures_filtered$signature_ct == "CL:0000840",]$signature
#
# x <- t(scores_mat)
# x <-  x[,sigs]
#
# # Remove rows wih NA
# x <- x[complete.cases(x), ]
# # Don't use fraction == 1
# x <- x[unlist(lapply(strsplit(rownames(x), "_"), "[", 2)) != "1",]
# # Check if CTOI not in controls
# controls <- unlist(lapply(strsplit(rownames(x), "_"), "[", 3))
# sum("CL:0000840" %in% controls)
#
# # Take only fractions of CTOI for y
# is_ctoi <- startsWith(rownames(x), "CL:0000840")
# y <- as.numeric(ifelse(is_ctoi, unlist(lapply(strsplit(rownames(x), "_"), "[", 2)), 0))
#
# # Lasso regression
# set.seed(123)
# res <- cv.glmnet(x, y, type.measure = "mse", alpha = 1)
# plot(res)
#
#
# p <- predict(res, newx = x, s = "lambda.min")
#
# lasso_coef <- res$glmnet.fit$beta[,res$glmnet.fit$lambda == res$lambda.min]
#
# cor(y, as.numeric(p))
#
# cor(y[1:25], p[1:25])
# cor(y[1:25], x[1:25,1])
#
# cor(y, p)
# cor(y, x[,11])
