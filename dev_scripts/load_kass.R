library(tidyverse)


celltype_conversion <- read_tsv("../xCell2.0/celltype_conversion_with_ontology.txt")
celltype_conversion_long <- celltype_conversion %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv", )
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene)
all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == gsub(pattern = "-", replacement = ".",all_models_annot$Sample))

# Tumor ref
tumor_ref <- all_models_expr[,!is.na(all_models_annot$Tumor_model_annot)]
tumor_labels <- all_models_annot[!is.na(all_models_annot$Tumor_model_annot),]

tumor_labels <- tumor_labels %>%
  left_join(celltype_conversion_long, by = c("Tumor_model_annot" = "all_labels")) %>%
  select(ont, label = xCell2_labels)

ontology_file_checked <- "../xCell2.0/kass_tumor_dependencies_checked.tsv"

ref <- as.matrix(tumor_ref)
labels <- as.data.frame(tumor_labels)
dim(ref)
dim(labels)

# Blood ref
blood_ref <- all_models_expr[,!is.na(all_models_annot$Blood_model_annot)]
blood_labels <- all_models_annot[!is.na(all_models_annot$Blood_model_annot),]

blood_labels <- blood_labels %>%
  left_join(celltype_conversion_long, by = c("Blood_model_annot" = "all_labels")) %>%
  select(ont, label = xCell2_labels)

ontology_file_checked <- "../xCell2.0/kass_blood_dependencies_checked.tsv"

ref <- as.matrix(tumor_ref)
labels <- as.data.frame(tumor_labels)
dim(ref)
dim(labels)
