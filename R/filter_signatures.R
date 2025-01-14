
filterSignatures <- function(ref, labels, pure_ct_mat, dep_list, signatures_collection, mixture_fractions, grubbs_cutoff, simulations_cutoff){

  # This function created mixtures for the simulations
  getMixtures <- function(ref, labels, ct, ct_data, dep_list, max_control_type,
                          mixture_fractions){

    getControlsMeanExpression <- function(ref, labels, ct, dep_list, max_control_type, mixture_fractions){

      controls2use <- names(dep_list)[!names(dep_list) %in% c(ct, dep_list[[ct]])]

      n_ds <- length(unique(labels$dataset)) # This number of datasets will effect the sampling of controls
      n_slice <- ifelse(n_ds >= max_control_type, 1, round((max_control_type/n_ds)+0.5))


      controls <- sapply(1:length(mixture_fractions), function(x){

        labels %>%
          filter(label %in% controls2use) %>% # Only controls
          group_by(dataset, label) %>%
          slice_sample(n=1) %>% # One cell type per dataset
          group_by(dataset) %>%
          slice_sample(n=n_slice) %>%
          ungroup() %>%
          slice_sample(n=max_control_type) %>%
          pull(sample) %>%
          ref[,.] %>%
          as.matrix() %>%
          Rfast::rowmeans()

        # labels %>%
        #   filter(label %in% controls2use) %>%
        #   group_by(dataset) %>%
        #   slice_sample(n=1) %>% # One cell type per dataset
        #   group_by(label) %>%
        #   slice_sample(n=1) %>% # One sample per datasets
        #   ungroup() %>%
        #   slice_sample(n=max_control_type) %>%
        #   pull(sample) %>%
        #   ref[,.] %>%
        #   as.matrix() %>%
        #   Rfast::rowmeans()

      })

      controls_fracs <- controls %*% diag(1-mixture_fractions)
      return(controls_fracs)

    }

    mixSmaples <- function(ref, labels, ct, dep_list, max_control_type, ct_mean_expression, mixture_fractions){

      # Generate a matrix of CTOI fractions:
      m <- matrix(rep(ct_mean_expression, length(mixture_fractions)), byrow = FALSE, ncol = length(mixture_fractions)) %*% diag(mixture_fractions)
      # Get random controls fractions
      c <- getControlsMeanExpression(ref, labels, ct, dep_list, max_control_type, mixture_fractions)
      # Add controls
      m <- m + c

      rownames(m) <- rownames(ref)
      colnames(m) <- as.character(mixture_fractions)

      return(m)

    }

    print(ct) #remove

    ct_data %>%
      group_by(dataset) %>%
      summarise(samples = list(sample)) %>%
      rowwise() %>%
      mutate(n_samples = length(samples)) %>%
      ungroup() %>%
      # Use top 20 references with most samples
      top_n(n = 20, wt = n_samples) %>%
      rowwise() %>%
      # Get cell type mean expression for each dataset
      mutate(ct_mean_expression = list(Rfast::rowmeans(as.matrix(ref[,samples])))) %>%
      mutate(mixtures = list(mixSmaples(ref, labels, ct, dep_list, max_control_type, ct_mean_expression, mixture_fractions))) %>%
      pull(mixtures) %>%
      return()
  }


  getMixturesCors <- function(signatures_collection, ct, mixture_ranked, mixture_fractions){

    sigs <- signatures_collection[startsWith(names(signatures_collection), paste0(ct, "#"))]

    # Score every ranked mixtures of CTOI
    cors <- sapply(mixture_ranked, function(ranked_mix){

      # Score
      scores <- sapply(sigs, simplify = TRUE, function(sig){
        singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
      })
      if (is.list(scores)) {
        sigs <- sigs[-which(lengths(scores) == 0)]
        scores <- sapply(sigs, simplify = TRUE, function(sig){
          singscore::simpleScore(ranked_mix, upSet = sig, centerScore = FALSE)$TotalScore
        })
      }
      colnames(scores) <- names(sigs)

      # Correlation
      apply(scores, 2, function(sig_scores){
        cor(sig_scores, mixture_fractions, method = "spearman")
      })

    })

    median_cors <- Rfast::rowMedians(cors)
    names(median_cors) <- rownames(cors)

    return(median_cors)
  }



  celltypes <- colnames(pure_ct_mat)

  scores_mat <- matrix(nrow = length(signatures_collection),
                       ncol = ncol(pure_ct_mat),
                       dimnames = list(names(signatures_collection), colnames(pure_ct_mat)))

  sig_type <- unlist(lapply(strsplit(names(signatures_collection), "#"), "[", 1))

  for (type in celltypes) {

    if (sum(type == sig_type) == 0) {
      errorCondition(paste0("No signatures found for cell type: ", type))
    }

    type_signatures <- signatures_collection[type == sig_type]
    dep_cells <- dep_list[[type]]

    if (length(dep_cells) > 0) {
      types_to_use <- !colnames(scores_mat) %in% dep_cells
    }else{
      types_to_use <- rep(TRUE, ncol(scores_mat))
    }

    sub_mix <- pure_ct_mat[,types_to_use]
    sub_mix_ranked <- singscore::rankGenes(sub_mix)
    for (i in 1:length(type_signatures)) {
      sig <- type_signatures[i]
      scores_out <- singscore::simpleScore(sub_mix_ranked, upSet = sig[[1]], centerScore = FALSE)$TotalScore
      scores_mat[which(rownames(scores_mat) == names(sig)), colnames(sub_mix_ranked)] <- scores_out
    }

  }

  scores_mat_tidy <- scores_mat %>%
    as_tibble(., rownames = NA) %>%
    rownames_to_column(var = "signature") %>%
    pivot_longer(cols = -signature, values_to = "score", names_to = "sample_ct") %>%
    separate(signature, into = "signature_ct", sep = "#", remove = FALSE, extra = "drop")%>%
    drop_na()

  # Filter by simulations
  simulations <- labels %>%
    group_by(ont, label) %>%
    nest() %>%
    rowwise() %>%
    mutate(mixture = list(getMixtures(ref = ref, labels = labels, ct = label, ct_data = data, dep_list = dep_list, max_control_type = 5, mixture_fractions = mixture_fractions))) %>%
    mutate(mixture_ranked = list(lapply(mixture, function(mix){singscore::rankGenes(mix)})))

  simulations.cors <- simulations %>%
    mutate(mixture_cor = list(getMixturesCors(signatures_collection = signatures_collection, ct = label, mixture_ranked = mixture_ranked, mixture_fractions = mixture_fractions)))

  simulations.filtered <- simulations.cors %>%
    mutate(sig_filtered = list(names(mixture_cor)[mixture_cor >= quantile(mixture_cor, simulations_cutoff)])) %>%
    pull(sig_filtered) %>%
    unlist()


  # Filter by Grubb's test
  #grubbs <- scores_mat_tidy %>%
  #  group_by(signature_ct, signature) %>%
  #  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  #  filter(grubbs_statistic >= quantile(grubbs_statistic, grubbs_cutoff)) %>%
  #  pull(signature)



  signatures_collection_filtered <- signatures_collection[names(signatures_collection) %in% simulations.filtered]

  filter_signature_out <- list("scoreMatTidy" = scores_mat_tidy, "sigCollectionFilt" = signatures_collection_filtered)

  return(filter_signature_out)
}
