library(tidyverse)
source("R/xCell2.R")

celltype_conversion_long <- read_tsv("Data/celltype_conversion_with_ontology.txt") %>%
  rowwise() %>%
  mutate(all_labels = str_split(all_labels, ";")) %>%
  unnest(cols = c(all_labels))

head()

######################## Single-cell RNA-seq references ---------------------------------------

############# Human (Tabula Sapiens)---------------------------------------


ts <- readRDS("/bigdata/almogangel/TabulaSapiens/tabula_sapiens.rds")
ts@assays$RNA@key <- "rna_"

# Use only 10X data
ts <- subset(x = ts, subset = assay == "10x 3' v3")

ts_labels <- tibble(ont = ts@meta.data$cell_type_ontology_term_id, label = ts@meta.data$cell_type,
                    sample = rownames(ts@meta.data), dataset = ts@meta.data$donor)
ts_ref <- ts@assays$RNA@counts


ref <- ts_ref
labels <- as.data.frame(ts_labels)
labels$ont <- as.character(labels$ont)
labels$label <- as.character(labels$label)

xCell2GetLineage(labels = labels[,1:2], out_file = "Data/ts_human_dependencies.tsv")

View(read.table("Data/ts_human_dependencies.tsv", sep = "\t", header = TRUE))
ontology_file_checked <- "Data/ts_human_dependencies.tsv"


########################  Bulk RNA-seq sorted cells references ---------------------------------------
############# Human ---------------------------------------
# Tumor ref (Kassandra) ----
# all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv")
all_models_expr <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
# all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == all_models_annot$Sample)


tumor_ref <- all_models_expr[,!is.na(all_models_annot$Tumor_model_annot)]
tumor_labels <- all_models_annot[!is.na(all_models_annot$Tumor_model_annot),]

tumor_labels <- tumor_labels %>%
  mutate(label = plyr::mapvalues(Tumor_model_annot, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)) %>%
  mutate(ont = plyr::mapvalues(label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)) %>%
  dplyr::select(ont, label, sample = Sample, dataset = Dataset)


ontology_file_checked <- "Data/kass_tumor_dependencies_checked.tsv"

ref <- as.matrix(tumor_ref)
labels <- as.data.frame(tumor_labels)
dim(ref)
dim(labels)

# Blood ref (Kassandra) ----
# all_models_expr <- read_csv("../xCell2.0/Kassandara_data/all_models_expr.tsv")
all_models_expr <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_expr.tsv")
all_models_expr <- data.frame(all_models_expr[,-1], row.names = all_models_expr$Gene, check.names = F)
# all_models_annot <- read_csv("../xCell2.0/Kassandara_data/all_models_annot.tsv")
all_models_annot <- read_csv("/bigdata/almogangel/kassandra_data/sorted_cells/all_models_annot.tsv")
colnames(all_models_annot)[1] <- "Sample"
all(colnames(all_models_expr) == all_models_annot$Sample)

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
laboratory_data_annotation <- read_tsv("/bigdata/almogangel/kassandra_data/sorted_cells/laboratory_data_annotation.tsv")
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
# saveRDS(blood_labels, "Data/kass_blood_labels_with_lab.rds")

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




# Super reference
sref <- readRDS("/bigdata/almogangel/super_ref_for_xcell2/human_bulk_ref.rds")

# Use main labels
sref_labels_main <- sref@labels[,c(1,2,5,6)]
colnames(sref_labels_main)[1:2] <- c("ont", "label")
xCell2GetLineage(sref_labels_main[,1:2], out_file = "/home/almogangel/xCell2_git/Data/sref_bulk_human_dependencies.tsv")

ontology_file_checked <- "Data/sref_bulk_human_dependencies_checked.tsv"

ref <- sref@data
labels <- sref_labels_main
dim(ref)
dim(labels)



######################## Microarry sorted cells references ---------------------------------------
############# Human ---------------------------------------
# LM222 -----------

lm22_ref <- read.table("../xCell2.0/LM22_source_GEPs.txt", header = TRUE, check.names = FALSE, row.names = 1, sep = "\t")

lm22_labels <- data.frame(sample = colnames(lm22_ref), dataset = rep("LM22", ncol(lm22_ref)))


lm22_labels$label <- plyr::mapvalues(lm22_labels$sample, celltype_conversion_long$all_labels, celltype_conversion_long$xCell2_labels, warn_missing = FALSE)
# T cells CD4 memory resting + T cells CD4 memory activated => CD4+ memory T-cells
lm22_labels[lm22_labels$label %in% c("T cells CD4 memory resting", "T cells CD4 memory activated"), ]$label <- "CD4+ memory T-cells"
# NK cells resting + NK cells activated = > NK cells
lm22_labels[lm22_labels$label %in% c("NK cells resting", "NK cells activated"), ]$label <- "NK cells"
# Dendritic cells resting + Dendritic cells activated => DC
lm22_labels[lm22_labels$label %in% c("Dendritic cells resting", "Dendritic cells activated"), ]$label <- "DC"
# Mast cells resting + Mast cells activated => Mast cells
lm22_labels[lm22_labels$label %in% c("Mast cells resting", "Mast cells activated"), ]$label <- "Mast cells"
lm22_labels$ont <- plyr::mapvalues(lm22_labels$label, celltype_conversion_long$xCell2_labels, celltype_conversion_long$ont, warn_missing = FALSE)


"B-cells" = lm22_labels[lm22_labels$label %in% c("naive B-cells", "Memory B-cells", "Plasma cells"),]
"CD4+ T-cells" = cbrx_lm22.out[c("CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells"),]),
"Macrophages" = cbrx_lm22.out[c("Macrophages M0", "Macrophages M1", "Macrophages M2"),]),
"T-cells" = cbrx_lm22.out[c("CD8+ T-cells", "CD4+ naive T-cells", "T cells CD4 memory resting", "CD4+ memory T-cells",
                                    "Follicular T-helper", "Tregs", "Tgd cells"),]), check.names = FALSE))


lm22_labels <- lm22_labels[,c(4,3,1,2)]

xCell2GetLineage(labels = lm22_labels[,1:2], out_file = "../xCell2.0/lm22_dependencies.tsv")
