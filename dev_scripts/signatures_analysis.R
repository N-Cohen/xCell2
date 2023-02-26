library(tidyverse)
library(gridExtra)

truths_dir <- "../xCell2.0/Kassandara_data/cell_values/"
mix_dir <- "../xCell2.0/Kassandara_data/expressions/"

ds <- gsub(".tsv", "", list.files(truths_dir))
# blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE115823", "GSE127813", "GSE53655", "GSE60424", "GSE64655", "sc_pbmc", "SDY67")
# tumor_ds <- ds[!ds %in% blood_ds]

corDF <- function(datasets = ds, ctoi, signatures_collection){

  signatures_ctoi <- signatures_collection[startsWith(names(signatures_collection), ctoi)]
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

    ctoi_frac <- truth[ctoi,]

    if (all(is.na(ctoi_frac))) {
      next
    }

    ctoi_frac <- ctoi_frac[1, !is.na(ctoi_frac)]
    ctoi_frac <- ctoi_frac[1, ctoi_frac != ""]

    if (all(ctoi_frac == 0)) {
      next
    }

    samples2use <- intersect(names(ctoi_frac), colnames(mix))
    ctoi_frac <- ctoi_frac[samples2use]
    mix <- mix[,samples2use]

    if (ncol(mix) < 4) {
      next
    }

    mix_ranked <- singscore::rankGenes(mix)


    if (!all(names(ctoi_frac) == colnames(mix_ranked))) {
      print(paste0("Problem with ds:" , file))
      break
    }


    # Score mixture
    scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
      singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
    })

    if (is.list(scores)) {
      signatures_ctoi <- signatures_ctoi[-which(lengths(scores) == 0)]
      scores <- sapply(signatures_ctoi, simplify = TRUE, function(sig){
        singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
      })
    }

    colnames(scores) <- names(signatures_ctoi)
    rownames(scores) <- colnames(mix_ranked)


    # Calculate correlations

    if (!all(names(ctoi_frac) == rownames(scores))) {
      print(paste0("Problem with ds:" , file))
      break
    }

    cors <- apply(scores, 2, function(x){
      cor(as.numeric(ctoi_frac), x, method = "spearman")
    })

    cors.df <- as.data.frame(cors)
    names(cors.df) <- file
    all_ds[[file]] <- cors.df
  }


  ctoi_cors.df <- bind_cols(all_ds)

  return(ctoi_cors.df)

}

getTopSigs <- function(ctoi, sig_col = signatures_collection_tumor){
  cors_sigs.df <- corDF(datasets = ds, ctoi = ctoi, signatures_collection = sig_col)
  top_sigs <- names(which(sort(apply(cors_sigs.df, 1, median), decreasing = TRUE) > 0.85))
  return(list(cors.df = cors_sigs.df, top_sigs = top_sigs))
}

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
