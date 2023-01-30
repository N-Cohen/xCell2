# xCell2Analysis function


xCell2Analysis <- function(counts, ref){

  # Score counts with ref's signatures
  counts_ranked <- singscore::rankGenes(counts)

  scores_out <- sapply(ref@signatures, function(x){
    singscore::simpleScore(counts_ranked, upSet=x, centerScore = FALSE)$TotalScore
  })
  colnames(scores_out) <- names(ref@signatures)
  rownames(scores_out) <- colnames(counts)

  # Make results tidy
  xCell2_out <- scores_out %>%
    as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to	= c("signature_ct", "quantile", "diff", "n_sim", "n_passed"), values_to = "score", names_sep = "_") %>%
    group_by(Sample, signature_ct) %>%
    summarise(mean_score = mean(score)) %>%
    left_join(ref@calibrationValues, by = "signature_ct") %>%
    mutate(adjusted_score = (mean_score-shift_value)*scale_factor) %>%
    select(Sample, signature_ct, adjusted_score) %>%
    pivot_wider(names_from = Sample, values_from = adjusted_score) %>%
    data.frame(row.names = 1)

  return(xCell2_out)

}





