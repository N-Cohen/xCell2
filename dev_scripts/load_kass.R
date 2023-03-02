library(tidyverse)

celltype_conversion <- read_tsv("Data/celltype_conversion_with_ontology.txt")
celltype_conversion_long <- celltype_conversion %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == all_models_annot$Sample)

# Tumor ref ----
tumor_ref <- all_models_expr[,!is.na(all_models_annot$Tumor_model_annot)]
tumor_labels <- all_models_annot[!is.na(all_models_annot$Tumor_model_annot),]

tumor_labels <- tumor_labels %>%
  mutate(label = plyr::mapvalues(Tumor_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label)


ontology_file_checked <- "../xCell2.0/kass_tumor_dependencies_checked.tsv"

ref <- as.matrix(tumor_ref)
labels <- as.data.frame(tumor_labels)
dim(ref)
dim(labels)

# Blood ref  ----
blood_ref <- all_models_expr[,!is.na(all_models_annot$Blood_model_annot)]
blood_labels <- all_models_annot[!is.na(all_models_annot$Blood_model_annot),]
all(colnames(blood_ref) == blood_labels[!is.na(blood_labels$Blood_model_annot),]$Sample)

blood_labels <- blood_labels %>%
  mutate(label = plyr::mapvalues(Blood_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label)

# Add lab data
laboratory_data_expressions <- read_tsv("../xCell2.0/Kassandara_data/laboratory_data_expressions.tsv")
laboratory_data_expressions <- data.frame(laboratory_data_expressions[,-1], row.names = laboratory_data_expressions$Gene, check.names = F)
laboratory_data_annotation <- read_tsv("../xCell2.0/Kassandara_data/laboratory_data_annotation.tsv")
colnames(laboratory_data_annotation)[1] <- "Sample"
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)
laboratory_data_expressions <- laboratory_data_expressions[,laboratory_data_annotation$Sample]
all(colnames(laboratory_data_expressions) == laboratory_data_annotation$Sample)

lab_labels <- laboratory_data_annotation %>%
  mutate(label = plyr::mapvalues(Cell_type, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label)

blood_labels <- rbind(blood_labels, lab_labels)
blood_ref <- cbind(blood_ref, laboratory_data_expressions)


ontology_file_checked <- "../xCell2.0/kass_blood_dependencies_checked.tsv"

ref <- as.matrix(blood_ref)
labels <- as.data.frame(blood_labels)
dim(ref)
dim(labels)




# BlueprintEncode ----


bp <- celldex::BlueprintEncodeData()

bp_ref <- as.matrix(bp@assays@data$logcounts)
bp_labels <- bp@colData

table(bp_labels$label.fine)
table(bp_labels$label.main)
bp_labels[bp_labels$label.main == "B-cells",]$label.fine <- "B-cells"


all(colnames(bp_ref) == rownames(bp_labels))

bp_labels <- bp_labels %>%
  as_tibble() %>%
  dplyr::select(ont = label.ont, label = label.fine) %>%
  mutate(label = plyr::mapvalues(label, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE))


sort(table(bp_labels$label)) # Minimum 3 samples for cell type
ct2keep <- names(which(table(bp_labels$label) > 2))
bp_ref <- bp_ref[,bp_labels$label %in% ct2keep]
bp_labels <- bp_labels[bp_labels$label %in% ct2keep,]

dim(bp_ref)
dim(bp_labels)

ontology_file_checked <- "../xCell2.0/bp_dependencies_checked.tsv"

ref <- as.matrix(bp_ref)
labels <- as.data.frame(bp_labels)
