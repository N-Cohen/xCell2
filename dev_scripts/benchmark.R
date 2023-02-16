library(tidyverse)

# xcell2_blood_refsigs <- xCell2Ref.s4@signatures
# xcell2_tumor_refsigs <- signatures_collection_filtered
# xcell2_bp_refsigs <- signatures_collection_filtered
blood_ds <- c("BG_blood", "GSE107011", "GSE107572", "GSE115823", "GSE127813", "GSE53655", "GSE60424", "GSE64655", "sc_pbmc", "SDY67")
#blood_ref <- c("EPIC_BRef", "LM22_matrix_with_BC_CIBERSORTX.tsv", "Scaden_PBMC", "xCell2.0 - Blood")
#blood_tumor <- c("EPIC_BRef", "LM22_matrix_with_BC_CIBERSORTX.tsv", "Scaden_PBMC", "xCell2.0 - Blood")


xCell.data <- xCell::xCell.data

truths_dir <- "../xCell2.0/Kassandara_data/cell_values/"
mix_dir <- "../xCell2.0/Kassandara_data/expressions/"
ds <- gsub(".tsv", "", list.files(truths_dir))

score_xCell2 <- function(mix, sigs){
  mix_ranked <- singscore::rankGenes(mix)
  scores <- sapply(sigs, simplify = TRUE, function(sig){
    singscore::simpleScore(mix_ranked, upSet = sig, centerScore = FALSE)$TotalScore
  })
  colnames(scores) <- names(sigs)
  rownames(scores) <- colnames(mix_ranked)
  scores <- as.data.frame(scores)
  scores <- scores %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(-sample, names_to = "signature", values_to = "score") %>%
    rowwise() %>%
    separate(signature, into = "celltype", sep = "#", remove = FALSE, extra = "drop") %>%
    group_by(sample, celltype) %>%
    summarise(score = mean(score)) %>%
    pivot_wider(names_from = sample, values_from = score) %>%
    as_data_frame()
  xcell2.out <- data.frame(scores[,-1], row.names = scores$celltype, check.names = FALSE)

}



getCorrelations <- function(dataset, celltypes2add = c("Monocytes", "Neutrophils", "NK cells", "cDC", "T-cells", "Tregs",
                                                       "Memory T-helpers", "CD4+ T-cells", "CD8+ T-cells", "B-cells",
                                                       "Fibroblasts", "Lymphocytes", "Cancer cells")){

  # Load mixture
  mix <- read.table(paste0("../xCell2.0/Kassandara_data/expressions/", dataset, "_expressions.tsv"), check.names = FALSE, row.names = 1, header = TRUE)

  # Load truth
  truth <- read_tsv(paste0("../xCell2.0/Kassandara_data/cell_values/", dataset, ".tsv")) %>%
    dplyr::rename(celltype = 1) %>% # Fix for some datasets
    filter(!endsWith(celltype, "l)")) %>%
    filter(celltype != "Respiratory_cells") %>%
    filter(celltype != "Tumor KI67+") %>%
    mutate(celltype = plyr::mapvalues(celltype, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE))

  samples2use <- intersect(colnames(truth)[-1], colnames(mix))
  mix <- mix[,samples2use]
  truth <- truth[,c(colnames(truth)[1], samples2use)]

  if(!all(colnames(mix) == colnames(truth)[-1])){
    print(paste0("Problem with ds:" , file))
    break
  }

  truth <- pivot_longer(truth, -celltype, names_to = "sample", values_to = "true_fracs")

  # Run xCell2
  xcell2_blood.out <- score_xCell2(mix, sigs = xcell2_blood_refsigs)
  xcell2_tumor.out <- score_xCell2(mix, sigs = xcell2_tumor_refsigs)
  if (dataset == "SC_ovarian_cancer") {
    xcell2_bp.out <- NULL
  }else{
    xcell2_bp.out <- score_xCell2(mix, sigs = xcell2_bp_refsigs)
  }

  # Run xCell
  xcell.out <- as.data.frame(xCell::xCellAnalysis(mix,rnaseq = FALSE))
  rownames(xcell.out) <- plyr::mapvalues(rownames(xcell.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  # Load results from other methods
  kass.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/Kassandra/", dataset ,"_predicted_by_Kassandra.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(kass.out) <- plyr::mapvalues(rownames(kass.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  abis.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/ABIS/", dataset, "_predicted_by_ABIS.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(abis.out) <- plyr::mapvalues(rownames(abis.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_bref.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/EPIC_BRef/", dataset, "_predicted_by_EPIC_BRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_bref.out) <- plyr::mapvalues(rownames(epic_bref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  epic_tref.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/EPIC_TRef/", dataset, "_predicted_by_EPIC_TRef.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(epic_tref.out) <- plyr::mapvalues(rownames(epic_tref.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  fardeep_abs.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/fardeep_absolute/", dataset, "_predicted_by_fardeep_absolute.tsv")
  if(file.exists(fardeep_abs.in)){
    fardeep_abs.out <- read.table(fardeep_abs.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_abs.out) <- plyr::mapvalues(rownames(fardeep_abs.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_abs.out <- NULL
  }
  fardeep_rel.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/fardeep_relative/", dataset, "_predicted_by_fardeep_relative.tsv")
  if(file.exists(fardeep_rel.in)){
    fardeep_rel.out <- read.table(fardeep_rel.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(fardeep_rel.out) <- plyr::mapvalues(rownames(fardeep_rel.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    fardeep_rel.out <- NULL
  }
  cbrx_hnsc.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/HNSC_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_HNSC_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(cbrx_hnsc.out) <- plyr::mapvalues(rownames(cbrx_hnsc.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  cbrx_lm22.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/LM22_matrix_with_BC_CIBERSORTX/", dataset, "_predicted_by_LM22_matrix_with_BC_CIBERSORTX.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  cbrx_lm22.out <- cbrx_lm22.out[!rownames(cbrx_lm22.out) %in% c("Mast cells resting", "Mast cells activated", "NK cells resting", "NK cells activated",
                                                                 "Dendritic cells resting", "Dendritic cells activated"),] # Fix for some datasets
  rownames(cbrx_lm22.out) <- plyr::mapvalues(rownames(cbrx_lm22.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/quantiseq/", dataset, "_predicted_by_quantiseq.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti.out) <- plyr::mapvalues(rownames(quanti.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  quanti_tumor.out <- read.table(paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/quantiseq_tumor/", dataset, "_predicted_by_quantiseq_tumor.tsv"), header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
  rownames(quanti_tumor.out) <- plyr::mapvalues(rownames(quanti_tumor.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)

  scaden.in <- paste0("../xCell2.0/Kassandara_data/predicted_by_algorithms/Scaden_PBMC/", dataset, "_predicted_by_Scaden_PBMC.tsv")
  if(file.exists(scaden.in)){
    scaden.out <- read.table(scaden.in, header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)
    rownames(scaden.out) <- plyr::mapvalues(rownames(scaden.out), celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
  }else{
    scaden.out <- NULL
  }


  all_methods <- c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode", "xCell", "Kassandara", "ABIS", "EPIC-BRef", "EPIC-TRef", "FARDEEP-Abs", "FARDEEP-Rel",
                   "CIBERSORTx-HNSC", "CIBERSORTx-LM22", "quanTIseq", "quanTIseq-T", "Scaden")

  all.out <- list(xcell2_blood.out, xcell2_tumor.out, xcell2_bp.out, xcell.out, kass.out, abis.out, epic_bref.out, epic_tref.out, fardeep_abs.out,
                  fardeep_rel.out, cbrx_hnsc.out, cbrx_lm22.out, quanti.out, quanti_tumor.out, scaden.out)
  names(all.out) <- all_methods

  all.out <- compact(all.out)
  shared_samples <- Reduce(intersect, lapply(all.out, names))
  shared_celltypes <- Reduce(intersect, lapply(all.out, rownames))
  shared_celltypes <- unique(c(shared_celltypes, celltypes2add))

  all.out <- lapply(all.out, function(x){x[shared_celltypes, shared_samples]})
  all.out <- lapply(all.out, function(x){cbind(celltype = shared_celltypes, x)})
  all.out <- lapply(all.out, function(x){drop_na(x)})
  all.out <- lapply(all.out, function(x){x[rowSums(x != 0) > 3,]}) # Minimum 3 values for cell type

  cors.out <- enframe(all.out, name = "method", value = "predictions") %>%
    unnest(cols = c(predictions)) %>%
    pivot_longer(-c(method, celltype), names_to = "sample", values_to = "prediction") %>%
    left_join(truth, by = c("celltype", "sample")) %>%
    mutate(true_fracs = as.numeric(true_fracs)) %>%
    drop_na() %>%
    group_by(method, celltype) %>%
    summarise(Pearson = cor(prediction, true_fracs),
              Spearman = cor(prediction, true_fracs, method = "spearman"),
              MSE = mean((prediction - true_fracs)^2)) %>%
    pivot_longer(cols = c("Pearson", "Spearman"), names_to = "cor_method", values_to = "cor_score") %>%
    mutate(dataset = dataset)

  return(cors.out)

}




cors.out.list <- vector("list", length(ds))
names(cors.out.list) <- ds

for (d in ds) {

#    if (d %in% ds[1:19]) {
#    next
#  }

  print(paste0(d, " --------------------------------------------------------- "))
  cors.out <- getCorrelations(dataset = d) %>%
    mutate(dataset = d)
  cors.out.list[[d]] <- cors.out
}


all_cors.out <- bind_rows(cors.out.list)

all_cors.out <- all_cors.out %>%
  rowwise() %>%
  mutate(ds_type = ifelse(dataset %in% blood_ds, "Blood datasets", "Tumor datasets"))


ctoi <- "Cancer cells"

all_cors.out %>%
  mutate(is_xcell2 = ifelse(method %in% c("xCell2.0 - Blood", "xCell2.0 - Tumor", "xCell2.0 - BlueprintEncode"), "yes", "no")) %>%
  ungroup() %>%
  filter(dataset != "Tonsils") %>%
  filter(cor_method == "Spearman") %>%
  filter(celltype == ctoi) %>%
  ggplot(., aes(x=method, y=cor_score, fill=is_xcell2)) +
  geom_boxplot(position = position_dodge(1), alpha = 0.6, outlier.shape = NA,  alpha = 0.5) +
  geom_point(aes(col=dataset), size=3, position = position_jitterdodge(jitter.width = .1, dodge.width = .5)) +
  facet_grid(. ~ ds_type) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(11, "Paired"),
                              "#424242", "#FFE7BA", "#8B1C62", "#CDC9C9", "#00F5FF", "#FF3E96")) +
  scale_fill_manual( values = c("yes"="tomato", "no"="gray"), guide = FALSE) +
  labs(y = "Correlation coefficient (r)", title = ctoi) +
  theme_linedraw() +
  theme(strip.background =element_rect(fill="#8B2323"),
        plot.title = element_text(size=22),
        axis.ticks = element_line(linetype = "blank"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "#5E5E5E",
                                        linetype = "dashed"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        strip.text.x = element_text(size = 16)) +
  labs(x = NULL, colour = NULL, fill = NULL)


cors.out.list[[1]] %>%
  filter(cor_method == "Spearman") %>%
  ggplot(., aes(x=method, y=cor_score)) +
  geom_boxplot(aes(fill=cor_method), position = position_dodge(1), alpha = 0.6, outlier.shape = NA) +
  geom_point(aes(col=celltype, group=cor_method), size=3,  position = position_jitterdodge(jitter.width = .3, dodge.width = 2)) +
  facet_grid(. ~ cor_method) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_fill_brewer(palette="Pastel1") +
  labs(y = "Correlation coefficient (r)") +
  theme(axis.ticks = element_line(linetype = "blank"),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_line(colour = "gray95",
                                        linetype = "dashed"), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 40, hjust=1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "gray97"),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        strip.text.x = element_text(size = 16)) +
  labs(x = NULL, colour = NULL, fill = NULL)
