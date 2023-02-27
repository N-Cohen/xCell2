library(tidyverse)

source("R/xCell2.R")
source("R/utils.R")

xCell.data <- xCell::xCell.data
run_xcell <- function(mixture, rnaseq){
  xCell_out <- xCell::xCellAnalysis(mixture, rnaseq = rnaseq)

  xCell_out <- xCell_out %>%
    as_tibble(rownames = "celltype") %>%
    pivot_longer(-celltype, names_to = "sample", values_to = "xCell")

  return(xCell_out)
}

run_xcell2 <- function(mix, ref){
  xCell2_out <- xCell2Analysis(mix, ref)

  xCell2_out <- xCell2_out %>%
    as_tibble(rownames = "celltype") %>%
    pivot_longer(-celltype, names_to = "sample", values_to = "xCell2") %>%
    mutate(celltype = gsub("-", "_", celltype))

  return(xCell2_out)
}

load_results <- function(file, results_by, is_truth = FALSE){
  res <- read_tsv(file)
  colnames(res)[1] <- "celltype"

  res.out <- res %>%
    pivot_longer(-celltype, names_to = "sample", values_to = ifelse(is_truth, "Truth", results_by))



  return(res.out)
}

merge_results <- function(results, truth, ct_to_use){

  # Merge all results
  results_merged <- truth %>%
    filter(celltype %in% ct_to_use) %>%
    mutate(Truth = as.numeric(Truth))

  for (res in results) {
    results_merged <- results_merged %>%
      left_join(res, by = c("celltype", "sample"))
  }

  # Calculate correlations
  results_merged <- results_merged %>%
    pivot_longer(., cols = -c(celltype, sample, Truth), names_to = "method", values_to = "score") %>%
    mutate(score = round(score, 3)) %>%
    #replace(is.na(.), 0) %>%
    drop_na() %>%
    group_by(celltype, method) %>%
    summarise("scores" = list(score), "true_fracs" = list(Truth)) %>%
    rowwise() %>%
    mutate(Pearson = cor(scores, true_fracs),
           Spearman = cor(scores, true_fracs, method = "spearman"),
           MSE = mean((true_fracs - scores)^2)) %>%
    mutate(Pearson = ifelse(is.na(Pearson), 0, Pearson),
           Spearman = ifelse(is.na(Spearman), 0, Spearman))

  return(results_merged)
}

plot_bars <- function(results.merged, ylimits = c(-1, 1)){

  results.merged %>%
    pivot_longer(cols = c("Pearson", "Spearman"), names_to = "cor_method", values_to = "cor_score") %>%
    ggplot(., aes(x=method, y=cor_score)) +
    geom_boxplot(aes(fill=cor_method), position = position_dodge(1), alpha = 0.6, outlier.shape = NA) +
    geom_point(aes(col=celltype, group=cor_method), size=3,  position = position_jitterdodge(jitter.width = .3, dodge.width = 2)) +
    facet_grid(. ~ cor_method) +
    scale_y_continuous(limits = ylimits) +
    scale_fill_brewer(palette="Pastel1") +
    labs(y = "Correlation coefficient (r)") +
    theme(axis.ticks = element_line(linetype = "blank"),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_line(colour = "gray95",
                                          linetype = "dashed"), axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          panel.background = element_rect(fill = "gray97"),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          strip.text.x = element_text(size = 16)) +
    labs(x = NULL, colour = NULL, fill = NULL)
}

plot_lines <- function(results.merged){

  methods <- unique(results.merged$method)
  p1 <- lapply(methods, function(x){
    p <- results.merged %>%
      filter(method == x) %>%
      unnest(cols = c(scores, true_fracs)) %>%
      ggplot(., aes(x=scores, y=true_fracs)) +
      geom_point() +
      stat_smooth(method = "lm") +
      facet_wrap(~celltype, scales = "free") +
      theme(axis.line = element_line(size = 0.2,
                                     linetype = "solid"), panel.grid.major = element_line(colour = NA),
            axis.text = element_text(size = 13),
            panel.background = element_rect(fill = NA),
            legend.background = element_rect(fill = NA)) +labs(x = "Scores", y = "True fraction", title = x)
    p
  })
  names(p1) <- methods

  gridExtra::grid.arrange(grobs = p1)

}

plot_spill <- function(results, truth, ct_to_use = NULL, cor_method = "spearman"){

  if (is.null(ct_to_use)) {
    ct_to_use <- unique(truth$celltype)
  }

  method_name <- colnames(results)[3]
  truth2scores <- truth %>%
    filter(celltype %in% ct_to_use) %>%
    left_join(results, by = c("celltype", "sample")) %>%
    drop_na()

  pairs <- crossing(ct_truth = unique(truth2scores$celltype), ct_score = unique(truth2scores$celltype)) %>%
    mutate(data = NA)

  for (i in 1:nrow(pairs)) {
    tmp <- truth2scores %>%
      filter(celltype == pull(pairs[i,1])) %>%
      select(sample, Truth) %>%
      mutate(Truth = as.numeric(Truth)) %>%
      full_join(select(filter(truth2scores, celltype == pull(pairs[i,2])), c(sample, method_name)), by = "sample") %>%
      drop_na()

    pairs[i,3] <- cor(pull(tmp[,2]), pull(tmp[,3]), method = cor_method)
  }

  pairs <- pivot_wider(pairs, names_from = ct_truth, values_from = data)
  M <- as.matrix(pairs[,-1])
  rownames(M) <- colnames(M)

  colors <- RColorBrewer::brewer.pal(11, "RdBu")
  corrplot::corrplot(M, method="circle", tl.cex = 1, tl.col = "black",
                     col=colorRampPalette(rev(colors))(200), order="original", title=paste0(method_name, " (", cor_method, ")"),
                     addCoef.col = "black", mar=c(0,0,1,0))

}

########################################################################################
# Build references
########################################################################################

# Blood ref ----
ct_ont <- read_tsv("/home/almogangel/xCell2/celltype_onthology.txt")

# Load data
# Blood model - Collection of 9056 bulk RNA-seq samples from 505 datasets of sorted cells, cancer cells and cell lines
all_models_expr <- read.table("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv", sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
# Load laboratory data - 348 bulk RNA-seq samples of sorted cell populations (including 343 samples of cells from the blood)
laboratory_data_expressions <- read.table("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_expressions.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
laboratory_data_annotation <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_annotation.tsv")
colnames(laboratory_data_annotation)[1] <- "Sample"


# Use blood-model cell types and match their ontology
labels_blood <- all_models_annot %>%
  dplyr::select(Sample, Blood_model_annot, Dataset) %>%
  drop_na() %>%
  mutate(ont = plyr::mapvalues(Blood_model_annot, ct_ont$label, ct_ont$ont)) %>%
  rename("label" = "Blood_model_annot")


# Use everything from the laboratory data
labels_lab <- laboratory_data_annotation %>%
  mutate(Dataset = "Lab") %>%
  mutate(ont = plyr::mapvalues(Cell_type, ct_ont$label, ct_ont$ont)) %>%
  rename("label" = Cell_type)


# There is only a single Basophils/Memory_CD8_T_cells sample in the entire data - remove
rm_baso_sample <- labels_blood$label != "Basophils"
labels_blood <- labels_blood[rm_baso_sample,]

rm_cd8_mem <- labels_lab$label != "Memory_CD8_T_cells"
labels_lab <- labels_lab[rm_cd8_mem,]
rm(rm_baso_sample, rm_cd8_mem)


# Some cell types are unique to the Lab data - split them for later
labels_lab_unique <- labels_lab[!labels_lab$label %in% labels_blood$label,]

# The non-unique Lab data can go to the blood model
labels_blood <- rbind(labels_blood, labels_lab[!labels_lab$Sample %in% labels_lab_unique$Sample,])


# Split blood-model data to train/test (66-33%)
celltypes_sorted <- names(sort(table(labels_blood$label), decreasing = TRUE))
ds_to_split <- c()
train_samples <- c()
for (ct in celltypes_sorted) {
  n_samples_for_train <- round(sum(labels_blood$label == ct)*0.66) # Number of samples to have ~66% samples of this cell-type
  ct_ds <- sort(table(labels_blood[labels_blood$label == ct,]$Dataset), decreasing = FALSE)
  ds_seen <- labels_blood[labels_blood$Sample %in% train_samples,]$Dataset   # Remove datasets that we have seen already
  ct_ds <- ct_ds[!names(ct_ds) %in% ds_seen]

  # If there is only one/none dataset save and skip
  if (length(ct_ds) == 1) {
    ds_to_split <- c(ds_to_split, names(ct_ds))
    next
  }else if(length(ct_ds) == 2){
    # Take the bigger dataset to train
    bigger_ds <- names(sort(ct_ds, decreasing = TRUE)[1])
    train_samples <- c(train_samples, labels_blood[labels_blood$Dataset == bigger_ds,]$Sample)
    next
  }

  for (i in 1:length(ct_ds)) {
    train_samples <- c(train_samples, labels_blood[labels_blood$Dataset == names(ct_ds[i]),]$Sample)  # Save all sample to train
    n_ct_in_train <- sum(labels_blood[labels_blood$Sample %in% train_samples,]$label == ct)  # Check how many cts are in train now
    if (n_ct_in_train >= n_samples_for_train) {
      break
    }
  }
}

# No ds_to_split
rm(ds_to_split)


# Split labels_lab_unique
unique_lab_cts <- unique(labels_lab_unique$label)
for (ct in unique_lab_cts) {
  ct_samples <- labels_lab_unique[labels_lab_unique$label == ct,]$Sample
  n_samples_for_train <- round(length(ct_samples)*0.66)
  train_samples <- c(train_samples, ct_samples[sample(1:length(ct_samples), n_samples_for_train)])
}


# Make sure train and test have all cell types
labels <- rbind(labels_blood, labels_lab_unique)
all_celltypes <- unique(labels$label)
all(all_celltypes %in% labels[labels$Sample %in% train_samples,]$label)
table(labels[labels$Sample %in% train_samples,]$label)
all(all_celltypes %in% labels[!labels$Sample %in% train_samples,]$label)

# Those cell types are missing from test - add them manually
missing_ct <- all_celltypes[!all_celltypes %in% labels[!labels$Sample %in% train_samples,]$label]

train_samples <- train_samples[!train_samples %in% labels[labels$Dataset == "GSE102215",]$Sample] # Granulocytes
train_samples <- train_samples[!train_samples %in% labels[labels$Dataset == "GSE107011",]$Sample] # Non_classical_monocytes & Classical_monocytes
train_samples <- train_samples[!train_samples %in% labels[labels$Dataset == "GSE81443",]$Sample] # Plasma_B_cells
train_samples <- train_samples[!train_samples %in% labels[labels$Dataset == "GSE115736",]$Sample] # Eosinophils
train_samples <- train_samples[!train_samples %in% labels[labels$Dataset == "GSE110999",]$Sample] # Naive_B_cells

# x <- labels %>%
#   filter(Sample %in% train_samples) %>%
#   group_by(Dataset, label) %>%
#   summarise(n=n()) %>%
#   filter(Dataset %in% labels[labels$label == missing_ct[6],]$Dataset)
# View(x)

# Check again
all(all_celltypes %in% labels[labels$Sample %in% train_samples,]$label)
all(all_celltypes %in% labels[!labels$Sample %in% train_samples,]$label)


# Ref data
ref <- cbind(all_models_expr, laboratory_data_expressions)
ref <- ref[,labels$Sample]
all(colnames(ref) == labels$Sample)

# Mark test samples
labels <- labels %>%
  mutate(is_test = ifelse(labels$Sample %in% train_samples, FALSE, TRUE)) %>%
  select(ont, label, is_test)


ref <- as.matrix(ref)
labels <- as.data.frame(labels)



# Tumor ref ------------
ct_ont <- read_tsv("Data/celltype_onthology.txt")
all_models_expr <- read.table("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv", sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"

labels <- all_models_annot %>%
  dplyr::select(Sample, Tumor_model_annot, Dataset) %>%
  drop_na() %>%
  mutate(ont = plyr::mapvalues(Tumor_model_annot, ct_ont$label, ct_ont$ont)) %>%
  rename("label" = Tumor_model_annot)

ref <- as.matrix(all_models_expr[,labels$Sample])

# Get cell-types ontology
ontology <- unique(as_tibble(labels[,c(2,4)]))
xCell2GetLineage(ontology_table = ontology, out_file = "Data/kass_tumor_dependencies.tsv")
ontology_file_checked <- "Data/kass_tumor_dependencies_checked.tsv"

# Split all_models_annot_tumor to train/test ~(66%-33%) - 8146 samples
train_ds <- c("E-MTAB-6643", "GSE100382", "GSE80727", names(sort(table(labels$Dataset), decreasing = T)[1:95]))
all_celltypes <- unique(labels$label)
all(all_celltypes %in% labels[labels$Dataset %in% train_ds,]$label)
all(all_celltypes %in% labels[!labels$Dataset %in% train_ds,]$label)

labels <- labels %>%
  mutate(is_test = ifelse(Dataset %in% train_ds, FALSE, TRUE)) %>%
  select(ont, label, is_test) %>%
  as.data.frame()

xcell2_kass_ref_tumor <- xCell2NewRef(ref, labels, data_type = "rnaseq")
saveRDS(xcell2_kass_ref_tumor, "/bigdata/almogangel/kassandra_data/kassandra_xcell2ref_tumor.rds")


########################################################################################
# Benchmark mixtures
########################################################################################


# --------- GSE107572 - quanTIseq (tumor) ---------
celltype_conversion <- read.table("/home/almogangel/xCell2_archive_jan23/celltype_conversion.txt", sep = "\t", header = TRUE, check.names = FALSE)
truth <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/GSE107572.tsv", is_truth = TRUE)
truth$celltype <- plyr::mapvalues(truth$celltype, from = celltype_conversion$GSE107572, to = celltype_conversion$label)
mix <- read.table("/bigdata/almogangel/kassandra_data/24_validation_datasets/expressions/GSE107572_expressions.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
mix <- as.matrix(mix)

# Run xCell
xcell.out <- run_xcell(mix, rnaseq = FALSE)
xcell.out$celltype <- plyr::mapvalues(xcell.out$celltype, from = celltype_conversion$xCell, to = celltype_conversion$label)

# Run xCell2
xcell2.out <- run_xcell2(mix, ref)

# Load results from other methods
kass.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/Kassandra/GSE107572_predicted_by_Kassandra.tsv", results_by = "Kassandra")
kass.out$celltype <- plyr::mapvalues(kass.out$celltype, from = celltype_conversion$Kassandra, to = celltype_conversion$label)
cbr.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/CIBERSORT/GSE107572_predicted_by_CIBERSORT.tsv", results_by = "CIBERSORT")
cbr.out$celltype <- plyr::mapvalues(cbr.out$celltype, from = celltype_conversion$CIBERSORT, to = celltype_conversion$label)
epic.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/EPIC_BRef/GSE107572_predicted_by_EPIC_BRef.tsv", results_by = "EPIC")
epic.out$celltype <- plyr::mapvalues(epic.out$celltype, from = celltype_conversion$EPIC, to = celltype_conversion$label)
quanti.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/quantiseq/GSE107572_predicted_by_quantiseq.tsv", results_by = "quanTIseq")
quanti.out$celltype <- plyr::mapvalues(quanti.out$celltype, from = celltype_conversion$quanTIseq, to = celltype_conversion$label)
abis.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/ABIS/GSE107572_predicted_by_ABIS.tsv", results_by = "ABIS")
abis.out$celltype <- plyr::mapvalues(abis.out$celltype, from = celltype_conversion$ABIS, to = celltype_conversion$label)

# Summarize results
results <- list(xcell.out, xcell2.out, kass.out, cbr.out, epic.out, quanti.out, abis.out)
unique(truth$celltype)
ct_to_use <- Reduce(intersect, lapply(results, function(x) pull(unique(x[,1]))))
results.merged <- merge_results(results, truth, ct_to_use)

# Plots
plot_bars(results.merged, ylimits = c(-1, 1))
plot_lines(results.merged)
plot_spill(results = results[[2]], truth, ct_to_use, cor_method = "spearman")

# --------- BG blood ---------
celltype_conversion <- read.table("/home/almogangel/xCell2_archive_jan23/celltype_conversion.txt", sep = "\t", header = TRUE, check.names = FALSE)
mix <- read.table("/bigdata/almogangel/kassandra_data/BG_blood/BG_blood_expressions.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Run xCell
xcell.out <- run_xcell(mix, rnaseq = FALSE)
xcell.out$celltype <- plyr::mapvalues(xcell.out$celltype, from = celltype_conversion$xCell, to = celltype_conversion$label)

# Run xCell2
xcell2.out <- run_xcell2(mix, ref)

# Load results from other methods
truth <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/cell_values/BG_blood.tsv", is_truth = TRUE)
kass.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/Kassandra/BG_blood_predicted_by_Kassandra.tsv", results_by = "Kassandra")
kass.out$celltype <- plyr::mapvalues(kass.out$celltype, from = celltype_conversion$Kassandra, to = celltype_conversion$label)
cbr.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/CIBERSORT/BG_blood_predicted_by_CIBERSORT.tsv", results_by = "CIBERSORT")
cbr.out$celltype <- plyr::mapvalues(cbr.out$celltype, from = celltype_conversion$CIBERSORT, to = celltype_conversion$label)
epic.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/EPIC_BRef/BG_blood_predicted_by_EPIC_BRef.tsv", results_by = "EPIC")
epic.out$celltype <- plyr::mapvalues(epic.out$celltype, from = celltype_conversion$EPIC, to = celltype_conversion$label)
quanti.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/quantiseq/BG_blood_predicted_by_quantiseq.tsv", results_by = "quanTIseq")
quanti.out$celltype <- plyr::mapvalues(quanti.out$celltype, from = celltype_conversion$quanTIseq, to = celltype_conversion$label)
abis.out <- load_results("/bigdata/almogangel/kassandra_data/24_validation_datasets/predicted_by_algorithms/ABIS/BG_blood_predicted_by_ABIS.tsv", results_by = "ABIS")
abis.out$celltype <- plyr::mapvalues(abis.out$celltype, from = celltype_conversion$ABIS, to = celltype_conversion$label)

# Summarize
results <- list(xcell.out, xcell2.out, kass.out, cbr.out, epic.out, quanti.out, abis.out)
ct_to_use <- Reduce(intersect, lapply(results, function(x) pull(unique(x[,1]))))
results.merged <- merge_results(results, truth, ct_to_use)

# Plots
plot_bars(results.merged, ylimits = c(0, 1))
plot_lines(results.merged)
plot_spill(results[[2]], truth, ct_to_use)

