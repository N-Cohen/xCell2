setwd("/Users/noam/Desktop/Technion/bioinformatics project/xCell2")

xCell2GetLineage <- function(labels, out_file){
  
  cl <- ontoProc::getCellOnto()
  
  labels_uniq <- labels %>%
    as_tibble() %>%
    unique()
  
  labels_uniq %>%
    rowwise() %>%
    mutate(descendants = list(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)),
           ancestors = list(ontologyIndex::get_ancestors(cl, terms = ont))) %>%
    mutate(ancestors = list(ancestors[ancestors != ont])) %>%
    mutate(descendants = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% descendants, 2]), collapse = ", "),
           ancestors = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% ancestors, 2]), collapse = ", ")) %>%
    write_tsv(out_file)
  
  warning("It is recommended that you manually check the cell-type ontology file: ", out_file)
}



##########################################################################################
# Train signatures
##########################################################################################

# NOTE (!!!):
#  - No underscore in cell types labels
#  - Cell type onthology with colon in example: CL:0000545 (not CL_0000545)

# UPDATES:
# - in createSignatures I made cor_cells_cutoff to be the top 10% value of all correlation values
# - per_cells have been removed and is now 0.5 !!!
# - cor_mat is now in spearman (not pearson)

# DEPENDENCIES:
# tidyverse, pheatmap, singscore, GSEABase, outliers, ontoProc, ontologyIndex

# USAGE:
# ref -  ref matrix of the new reference (genes x samples) (!!!)
# (!!!) Check that the numbers are not strings.

# labels - a data frame with rows correspond to samples in ref:
#   (1) first column for cell type onthology
#   (2) second column for cell type name (charaters)


# Remove
if (1 == 0) {
  data_type = "rnaseq"; mixture_fractions = c(0.001, 0.005, seq(0.01, 0.25, 0.02))
  probs = c(.1, .25, .33333333, .5); diff_vals = c(0, 0.1, 0.585, 1, 1.585, 2, 3, 4, 5)
  min_genes = 5; max_genes = 500; is_10x = FALSE
}
