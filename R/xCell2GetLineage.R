setwd("/Users/noam/Desktop/Technion/bioinformatics project/xCell2")
#########################################################################################
# Generate cell-type lineage automatically for dependencies
#########################################################################################

# USAGE:
# (1) labels - a data frame with:
#   (a) first column for cell type onthology - named "ont" (character)
#   (b) second column for cell type name - named "label" (character)
# (2) out_file - path to cell types lineage file for manual check

# DEPENDENCIES:
# tidyverse, ontoProc, ontologyIndex

xCell2GetLineage <- function(labels, out_file){

  cl <- ontoProc::getCellOnto()

  labels_uniq <- labels %>%
    as_tibble() %>%
    select(ont, label) %>%
    unique()

  labels_uniq %>%
    rowwise() %>%
    mutate(descendants = list(ontologyIndex::get_descendants(cl, roots = ont, exclude_roots = TRUE)),
           ancestors = list(ontologyIndex::get_ancestors(cl, terms = ont))) %>%
    mutate(ancestors = list(ancestors[ancestors != ont])) %>%
    mutate(descendants = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% descendants, 2]), collapse = ";"),
           ancestors = paste(pull(labels_uniq[pull(labels_uniq[,1]) %in% ancestors, 2]), collapse = ";")) %>%
    write_tsv(out_file)

  warning("It is recommended that you manually check the cell-type ontology file: ", out_file)
}
