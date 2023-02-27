library(tidyverse)

# Build a Signature Matrix File for CIBERSORTx
tumor_ref
tumor_labels

colnames(tumor_ref) <- tumor_labels$label
tumor_ref %>%
  rownames_to_column(var = "genes") %>%
  write_tsv("../xCell2.0/cibersortx_tumor_sig_mat.txt")

all_celltypes <- unique(tumor_labels$label)
pheno.df <- matrix(rep(1, ncol(tumor_ref)*length(all_celltypes)), ncol = ncol(tumor_ref))
pheno.df <- cbind(unique(tumor_labels$label), pheno.df)

pheno.df <- t(sapply(all_celltypes, function(ct){

  dep_vec <- sapply(tumor_labels$label, function(label){
    ifelse(ct == label, 1,
           ifelse(label %in% dep_list[[ct]], 0, 2))
  })

  as.numeric(dep_vec)

}))

pheno.df %>%
  as.data.frame(.) %>%
  rownames_to_column() %>%
  write_tsv("../xCell2.0/cibersortx_tumor_pheno_df.txt", col_names = FALSE)



# Run CIBERSORTx
run_cibersortx <- function(mix, sig_mat, dir, method="fractions", single_cell=FALSE, QN=FALSE,
                           absolute=FALSE, rmbatchBmode=FALSE, rmbatchSmode=FALSE){

  # Move mixture and signature matrix to CIBERSORTx directory
  file.copy(from = mix, to = paste0(dir, "/mixture.txt"), overwrite = TRUE)
  file.copy(from = sig_mat, to = paste0(dir, "/sigmat.txt"), overwrite = TRUE)

  # Make results directory
  results_dir <- paste0(dir, "/results")

  # Note: Currently I only use the "fractions" method of CIBERSORTx
  token <- "5f0976eaee0e5a995dda7b09596d1d4e"

  # Clean old results
  if(length(list.files(results_dir)) > 0){
    system(paste0("rm ", results_dir, "/*"))
  }

  cmd <- paste0("cd ", dir, "; ",
                "docker run -v ", dir, ":/src/data -v ", results_dir, ":/src/outdir cibersortx/",
                method, " --username almog.angel@campus.technion.ac.il --token ", token, " --single_cell ", single_cell,
                " --sigmatrix sigmat.txt --mixture mixture.txt --QN ", QN, " --absolute ",
                absolute, " --rmbatchBmode ", rmbatchBmode, " --rmbatchSmode ", rmbatchSmode,
                " 1> results/cibersortx.stdout 2> results/cibersortx.stderr")

  # Run Docker via shell
  system(cmd, wait = TRUE)

  # Load results
  res_file <- ifelse(rmbatchBmode | rmbatchSmode, "/CIBERSORTx_Adjusted.txt", "/CIBERSORTx_Results.txt")
  res <- read.table(paste0(results_dir, res_file), stringsAsFactors = FALSE, sep='\t', header = TRUE, row.names=1, check.names = FALSE)
  res <- t(res[,1:(ncol(res)-3)])

  cibersortx_out <- res %>%
    as_tibble(rownames = "celltype") %>%
    pivot_longer(-celltype, names_to = "sample", values_to = "CIBERSORTx")

  return(cibersortx_out)
}

makeCIBERSORTSigMat <- function(ref = "../xCell2.0/cibersortx_tumor_sig_mat.txt",
                                pheno = "../xCell2.0/cibersortx_tumor_pheno_df.txt",
                                sig_out = "../xCell2.0/cibersortx_tumor_ref.txt", method = "fractions", QN = FALSE, single_cell = FALSE){


  token <- "5f0976eaee0e5a995dda7b09596d1d4e"

  cmd <- paste0("cd ", dir, "; ",
                "docker run -v ", dir, ":/src/data -v ", dir, ":/src/outdir cibersortx/",
                method, " --username almog.angel@campus.technion.ac.il --token ", token,
                " --single_cell FALSE --refsample ref.txt --phenoclasses pheno.txt --QN ", QN, " 1> stdout 2> stderr")

  system(cmd, wait = TRUE)


  print("Done!")
}
