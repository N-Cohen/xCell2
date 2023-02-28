library(tidyverse)
library(gridExtra)

truths_dir <- "../xCell2.0/Kassandara_data/cell_values/"
mix_dir <- "../xCell2.0/Kassandara_data/expressions/"
validation_ds <- gsub(".tsv", "", list.files(truths_dir))[c(1:10, 12, 22)]
validation_ds_blood <- validation_ds[!validation_ds %in% c("ccRCC_cytof_CD45+", "NSCLC_cytof", "WU_ccRCC_RCCTC", "GSE120444", "GSE115823", "GSE121127")]


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
  rownames(scores) <- colnames(mix_ranked)
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
    truth <- truth[!endsWith(rownames(truth), "l)"),] # Fix for some datasets
    truth <- truth[rownames(truth) != "Respiratory_cells",] # Fix for some datasets
    truth <- truth[rownames(truth) != "Tumor KI67+",] # Fix for some datasets
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

  cos_final <- tibble(id = names(all_ds_cors), cor = all_ds_cors) %>%
    separate(id, into = c("dataset", "celltype", "signature"), sep = "\\.", extra = "merge") %>%
    rowwise() %>%
    mutate(passed_filter = ifelse(signature %in% filtered_sigs, "yes", "no"))


  return(cos_final)

}

getTopSigs <- function(ctoi, sig_col = signatures_collection_tumor){
  cors_sigs.df <- corDF(datasets = ds, ctoi = ctoi, signatures_collection = sig_col)
  top_sigs <- names(which(sort(apply(cors_sigs.df, 1, median), decreasing = TRUE) > 0.85))
  return(list(cors.df = cors_sigs.df, top_sigs = top_sigs))
}


# Get correlations for all blood ref signatures
xcell2_blood_ref <- xCell2Train(ref, labels, ontology_file_checked,
                                data_type = "rnaseq",  score_method = "singscore", mixture_fractions = c(0.001, seq(0.01, 0.25, 0.02), 1),
                                probs = c(.1, .25, .33333333, .5), diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5),
                                min_genes = 5, max_genes = 500)

blood_ref_sigs_cors <- getSigsCor(datasets = validation_ds_blood, signatures_collection = xcell2_blood_ref@all_signatures, filtered_sigs = names(xcell2_blood_ref@filtered_signatures))
# blood_ref_sigs_cors <- cos_final
blood_ref_sigs_cors %>%
  ungroup() %>%
  mutate(passed_filter = factor(blood_ref_sigs_cors$passed_filter, levels = c("yes", "no"))) %>%
  ggplot(., aes(x=dataset, y=cor)) +
  geom_violin(position = position_dodge(1), alpha = 0.6, fill="#C1CDCD") +
  geom_jitter(aes(col=passed_filter), size=2, alpha=0.6, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  scale_color_manual(values=c("#66CD00", "#CD3333")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "#008B8B", size=0.5) +
  facet_wrap(~celltype, scales = "free_x") +
  labs(y = "Spearman r", x = "")










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
