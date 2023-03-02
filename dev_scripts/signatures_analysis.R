library(tidyverse)
library(gridExtra)

truths_dir <- "../xCell2.0/validation_data//cell_values/"
mix_dir <- "../xCell2.0/validation_data/expressions/"
validation_ds <- gsub(".tsv", "", list.files(truths_dir))
validation_ds_blood <- c("BG_blood", "GSE107011", "GSE107572", "GSE127813", "GSE53655", "GSE60424", "sc_pbmc", "SDY67", "SDY420", "SDY311", "DREAM")


scoreMixtures <- function(ctoi, mixture_ranked, signatures_collection){
  signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), paste0(ctoi, "#"))]

  scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
    singscore::simpleScore(mixture_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })

  # In case some signatures contain genes that are all not in the mixtures
  if (is.list(scores)) {
    signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })
  }
  colnames(scores) <- names(signatures_ctoi)
  rownames(scores) <- colnames(mixture_ranked)
  return(t(scores))
}

getSigsCor <- function(datasets, signatures_collection, filtered_sigs){

  all_celltypes <- unique(gsub(pattern = "#.*", "", names(signatures_collection)))
  all_ds <- sapply(datasets, function(x) NULL)

  for (file in datasets) {

    # Load mixture
    mix <- read.table(paste0(mix_dir, file, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

    # Load truth
    truth <- read.table(paste0(truths_dir, file, ".tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(truth) <- plyr::mapvalues(rownames(truth), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

    samples2use <- intersect(colnames(truth), colnames(mix))
    mix <- mix[,samples2use]
    truth <- truth[,samples2use]

    mix_ranked <- singscore::rankGenes(mix)

    if (!all(colnames(truth) == colnames(mix_ranked))) {
      errorCondition(paste0("Error with dataset: ", file))
    }


    all_celltypes_cor <- sapply(all_celltypes, function(ctoi){

      scores_ctoi <- scoreMixtures(ctoi, mix_ranked, signatures_collection)
      truth_ctoi <- truth[ctoi, colnames(scores_ctoi)]
      truth_ctoi <- truth_ctoi[1, !is.na(truth_ctoi)]
      truth_ctoi <- truth_ctoi[1, truth_ctoi != ""]
      scores_ctoi <- scores_ctoi[,names(truth_ctoi)]

      if (!all(names(truth_ctoi) == colnames(scores_ctoi))) {
        errorCondition(paste0("Error with dataset: ", file))
      }

      truth_ctoi <- as.numeric(truth_ctoi)

      if (all(truth_ctoi == 0)) {
        NULL
      }else{
        apply(scores_ctoi, 1, function(x){cor(x, truth_ctoi, method = "spearman")})
      }


    })

    all_ds[[file]] <- all_celltypes_cor
  }

  all_ds_cors <- unlist(all_ds)

  cor_final <- tibble(id = names(all_ds_cors), cor = all_ds_cors) %>%
    separate(id, into = c("dataset", "celltype", "signature"), sep = "\\.", extra = "merge") %>%
    rowwise() %>%
    mutate(passed_filter = ifelse(signature %in% filtered_sigs, "yes", "no"))


  return(cor_final)

}

getTopSigs <- function(cor_final, top = 0.1){
  top_sigs <- blood_ref_sigs_cors %>%
    group_by(celltype, signature) %>%
    summarise(median_cor = median(cor)) %>%
    filter(median_cor >= quantile(median_cor, 1-top)) %>%
    pull(signature)

  return(top_sigs)
}


# Get correlations for all blood ref signatures
xcell2_blood_ref <- xCell2Train(ref, labels, ontology_file_checked,
                                data_type = "rnaseq",  score_method = "singscore", mixture_fractions = c(0.001, seq(0.01, 0.25, 0.02), 1),
                                probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                                min_genes = 5, max_genes = 500)

# saveRDS(xcell2_blood_ref, "../xCell2.0/tmp_R_data/blood_ref_1.3.23_grubb07.rds")

# Filter by Grubb's test - all
grubbs.filtered <- xcell2_blood_ref@score_mat %>%
  group_by(signature_ct, signature) %>%
  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  filter(grubbs_statistic >= quantile(grubbs_statistic, 0.7)) %>%
  pull(signature)

# Filter by Grubb's test - similarity

ct_similarity <- xcell2_blood_ref@score_mat %>%
  rowwise() %>%
  mutate(similarity = xcell2_blood_ref@correlationMatrix[signature_ct, sample_ct]) %>%
  select(signature_ct, sample_ct, similarity) %>%
  unique() %>%
  group_by(signature_ct) %>%
  mutate(similarity_level = ifelse(similarity >= quantile(similarity, 0.8), "high", "low")) %>%
  mutate(similarity_level = ifelse(signature_ct == sample_ct, "same", similarity_level))


grubbs.filtered.sim<- xcell2_blood_ref@score_mat %>%
  filter(signature %in% grubbs.filtered) %>%
  left_join(ct_similarity, by = c("signature_ct", "sample_ct")) %>%
  filter(similarity_level != "high") %>%
  group_by(signature_ct, signature) %>%
  summarise(grubbs_statistic = outliers::grubbs.test(score, type = 10, opposite = FALSE, two.sided = FALSE)$statistic[1]) %>%
  filter(grubbs_statistic >= quantile(grubbs_statistic, 0.9)) %>%
  pull(signature)




blood_ref_sigs_cors <- getSigsCor(datasets = validation_ds_blood, signatures_collection = xcell2_blood_ref@all_signatures, filtered_sigs = grubbs.filtered.sim)

grubbs
grubbs.high
grubbs.low <- blood_ref_sigs_cors

grubbs.low <- blood_ref_sigs_cors %>%
  group_by(dataset, celltype, passed_filter) %>%
  summarise(median_cor = median(cor)) %>%
  group_by(dataset, celltype) %>%
  arrange(passed_filter, .by_group = T) %>%
  mutate(delta_median =  median_cor - lag(median_cor)) %>%
  drop_na() %>%
  select(dataset, celltype, delta_median) %>%
  group_by(celltype) %>%
  summarise(delta_median = median(delta_median)) %>%
  arrange(-delta_median)


blood_ref_sigs_cors %>%
  ungroup() %>%
  filter(celltype %in% unique(blood_ref_sigs_cors$celltype)) %>%
  mutate(passed_filter = factor(passed_filter, levels = c("yes", "no"))) %>%
  ggplot(., aes(x=dataset, y=cor)) +
  geom_boxplot(aes(fill=passed_filter)) +
  scale_fill_manual(values=c("#66CD00", "#CD3333")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~celltype, scales = "free_x") +
  labs(y = "Spearman r", x = "") +
  theme(axis.text.x = element_text(angle = 40, hjust=1, face = "bold"))










  # All top sigs ----
cts <- c("T-cells", "CD4+ T-cells", "CD8+ T-cells", "B-cells",
         "Fibroblasts", "Cancer cells")
all_top.list <- lapply(cts, function(ct){getTopSigs(ct)})


# B-cells ----
b_cors_blood_sigs.df <- corDF(datasets = ds, ctoi = "B-cells", signatures_collection = signatures_collection_blood)
b_cors_tumor_sigs.df <- corDF(datasets = ds, ctoi = "B-cells", signatures_collection = signatures_collection_tumor)

p1 <- pheatmap::pheatmap(b_cors_blood_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]
p2 <- pheatmap::pheatmap(b_cors_tumor_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]

gl <- list("Blood Ref" = p1, "Tumor Ref" = p2)
grid.arrange(p1, p2, ncol=2, top="Blood Ref                                                                                                               Tumor Ref")


# Choose best sings - Blood
b_cors_blood_sigs.df_filtered <- b_cors_blood_sigs.df

top_sigs <- names(which(sort(apply(b_cors_blood_sigs.df_filtered, 1, median), decreasing = TRUE) > 0.85))

top_sigs <- names(sort(apply(b_cors_blood_sigs.df_filtered[top_sigs,], 1, mean), decreasing = TRUE))
pheatmap::pheatmap(b_cors_blood_sigs.df_filtered[top_sigs,], cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))


# CD8+ T-cells ----
cd8_cors_blood_sigs.df <- corDF(datasets = ds, ctoi = "CD8+ T-cells", signatures_collection = signatures_collection_blood)
cd8_cors_tumor_sigs.df <- corDF(datasets = ds, ctoi = "CD8+ T-cells", signatures_collection = signatures_collection_tumor)

p1 <- pheatmap::pheatmap(cd8_cors_blood_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]
p2 <- pheatmap::pheatmap(cd8_cors_tumor_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]

gl <- list("Blood Ref" = p1, "Tumor Ref" = p2)
grid.arrange(p1, p2, ncol=2, top="Blood Ref                                                                                                               Tumor Ref")


# Choose best sings - Blood
cd8_cors_blood_sigs.df_filtered <- cd8_cors_blood_sigs.df

top_sigs <- names(which(sort(apply(cd8_cors_blood_sigs.df_filtered, 1, median), decreasing = TRUE) > 0.85))

top_sigs <- names(sort(apply(cd8_cors_blood_sigs.df_filtered[top_sigs,], 1, mean), decreasing = TRUE))
pheatmap::pheatmap(cd8_cors_blood_sigs.df_filtered[top_sigs,], cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))

# CD4+ T-cells ----
cd4_cors_blood_sigs.df <- corDF(datasets = ds, ctoi = "CD4+ T-cells", signatures_collection = signatures_collection_blood)
signatures_collection_tumor <- signatures_collection
cd4_cors_tumor_sigs.df <- corDF(datasets = ds, ctoi = "CD4+ T-cells", signatures_collection = signatures_collection_tumor)

p1 <- pheatmap::pheatmap(cd4_cors_blood_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]
p2 <- pheatmap::pheatmap(cd4_cors_tumor_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))[[4]]

gl <- list("Blood Ref" = p1, "Tumor Ref" = p2)
grid.arrange(p1, p2, ncol=2, top="Blood Ref                                                                                                               Tumor Ref")


# Choose best sings - Blood
cd4_cors_blood_sigs.df_filtered <- cd4_cors_blood_sigs.df

top_sigs <- names(which(sort(apply(cd4_cors_blood_sigs.df_filtered, 1, median), decreasing = TRUE) > 0.85))

top_sigs <- names(sort(apply(cd4_cors_blood_sigs.df_filtered[top_sigs,], 1, mean), decreasing = TRUE))
pheatmap::pheatmap(cd4_cors_blood_sigs.df_filtered[top_sigs,], cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))

# Choose best sings - Tumor
cd4_cors_tumor_sigs.df_filtered <- cd4_cors_tumor_sigs.df

top_sigs <- names(which(sort(apply(cd4_cors_tumor_sigs.df_filtered, 1, median), decreasing = TRUE) > 0.85))

top_sigs <- names(sort(apply(cd4_cors_tumor_sigs.df_filtered[top_sigs,], 1, mean), decreasing = TRUE))
pheatmap::pheatmap(cd4_cors_tumor_sigs.df_filtered[top_sigs,], cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))


# CD8+ T-cells PD1 low ----

pd1_cors_tumor_sigs.df <- corDF(datasets = ds, ctoi = "CD8+ T-cells PD1 high", signatures_collection = signatures_collection_filtered)
pheatmap::pheatmap(pd1_cors_tumor_sigs.df, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                         breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))



# Manually ----
unique(tumor_labels$label)
xcell2_tumor_refsigs
x <- corDF(datasets = ds, ctoi = "CD8+ T-cells PD1 high", xcell2_tumor_refsigs)
pheatmap::pheatmap(x, cluster_rows=F, cluster_cols=F, scale = "none", col= c(RColorBrewer::brewer.pal(9, "YlOrRd"), "black"),
                   breaks = seq(0,1,0.1), legend_breaks  = seq(0,1,0.1))
