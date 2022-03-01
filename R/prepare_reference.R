## Prepare reference:
# Input: RNA-seq or Array count matrix with cell types in column and genes in rows
# Output: SummarizedExperiment class ready to use for xCell2

# (1) Cell types (column) - Write a script that try to convert cell types automatically.
#                           If he is not sure - he ask the user to enter annotation manually.
# (2) Counts -  Should we normalize all counts the same way?
# (3) Genes - Convert between gene formats if nessasary
# (4) Data type - Choose one: RNA-seq / Array
# (5) SummarizedExperiment - Store everything in a SE class



prepareRef <- function(ref){





}


## Take shared genes (for multiple references)
# Input: list of SummarizedExperiment class references
# Output: list of SummarizedExperiment class references with shared genes

sharedGenes <- function(ref_list){




}
