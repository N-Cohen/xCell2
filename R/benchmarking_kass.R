source("R/xCell2.R")
library(tidyverse)

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





# # Train (5364)
# sum(sort(table(all_models_annot_tumor$Dataset), decreasing = T)[1:95]) # 5364
# train_ds <- names(sort(table(all_models_annot_tumor$Dataset), decreasing = T)[1:95])
# all_models_annot_tumor_train <- all_models_annot_tumor[all_models_annot_tumor$Dataset %in% train_ds,]
#
# # Test (2782)
# all_ds <- unique(all_models_annot_tumor$Dataset)
# test_ds <- all_ds[!all_ds %in% train_ds]
# nrow(all_models_annot_tumor[all_models_annot_tumor$Dataset %in% test_ds,])
# all_models_annot_tumor_test <- all_models_annot_tumor[all_models_annot_tumor$Dataset %in% test_ds,]
#
# # Check if all cell-types exist in test and train
# table(all_models_annot_tumor_train$Tumor_model_annot)
# table(all_models_annot_tumor_test$Tumor_model_annot)
# # Macrophages_M2 are missing in train
# unique(all_models_annot_tumor_test$Tumor_model_annot)[!unique(all_models_annot_tumor_test$Tumor_model_annot) %in% unique(all_models_annot_tumor_train$Tumor_model_annot)]
# table(all_models_annot_tumor_test[all_models_annot_tumor_test$Tumor_model_annot == "Macrophages_M2",]$Dataset)
# # Move E-MTAB-6643, GSE100382 and GSE80727 to train
# ds_to_move <- c("E-MTAB-6643", "GSE100382", "GSE80727")
# all_models_annot_tumor_train <- rbind(all_models_annot_tumor_train,
#                                       all_models_annot_tumor_test[all_models_annot_tumor_test$Dataset %in% ds_to_move,])
# all_models_annot_tumor_test <- all_models_annot_tumor_test[!all_models_annot_tumor_test$Dataset %in% ds_to_move,]
#
#
# # Subset count matrix into test and train
# all_models_expr_tumor_train <- all_models_expr[,all_models_annot_tumor_train$Sample]
# all_models_expr_tumor_test <- all_models_expr[,all_models_annot_tumor_test$Sample]
# rm(all_models_expr)
#
# # Prepare train
# all(all_models_annot_tumor_train$Sample == colnames(all_models_expr_tumor_train))
# colnames(all_models_annot_tumor_train)[2] <- "Cell_type"
# ref <- all_models_expr_tumor_train
# labels <- all_models_annot_tumor_train
# all(labels$Sample == colnames(ref))
# labels <- labels %>%
#   dplyr::select(ont, Cell_type) %>%
#   dplyr::rename("label" = Cell_type)
# ref <- as.matrix(ref)
# labels <- as.data.frame(labels)
#
# # Prepare test
# all(all_models_annot_tumor_test$Sample == colnames(all_models_expr_tumor_test))
# colnames(all_models_annot_tumor_test)[2] <- "Cell_type"
# ref_test <- all_models_expr_tumor_test
# labels_test <- all_models_annot_tumor_test
# all(labels_test$Sample == colnames(ref_test))
# labels_test <- labels_test %>%
#   dplyr::select(ont, Cell_type) %>%
#   dplyr::rename("label" = Cell_type)
# ref_test <- as.matrix(ref_test)
# labels_test <- as.data.frame(labels_test)
#
