library(tidyverse)

# Script to create reference based on Mahmoud and Kassandra's bulk RNA-seq sorted cells data

sref <- readRDS("/bigdata/mahmoudy/humanCellTypeAtlas.rds")

dim(sref$data)

unique(sref$tissue_titles)

# Get main and fine tissues from "tissue_titles"
table(sref$lowest.label)
table(sref$tissue_titles)
tibble(tissue_main = NA, tissue_fine = NA, tissue_info = sref$tissue_titles, label = sref$lowest.label,
       is_cancer = sref$is_cancer, is_cell_line = sref$is_cell_line, ont = sref$cl.id, sample = as.character(sref$samples), dataset = as.character(sref$series)) %>%
  write_tsv("../xCell2/super_ref_tissues.tsv")

# Get labels tibble - <ont> <label> <sample> <dataset>
labels.tbl <- tibble(ont = sref$cl.id, label = sref$lowest.label, sample = sref$samples, dataset = sref$series)


sref.tbl <- tibble(organism = "human", is_cancer = sref$is_cancer, is_cell_line = sref$is_cell_line,
                   tissue_main = , tissue_fine = , tissue_info = sref$tissue_titles, labels =)
