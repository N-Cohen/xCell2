library(singscore)
library(SummarizedExperiment)
library(GSEABase)
library(ggplot2)
library(ggrepel)

get_label_list <- function(ref_list) {
  lapply(ref_list, function(ref) {
    levels(as.factor(ref$label.fine))
  })
}


reference <- readRDS(file = "~/Documents/xCell2.0/xCell2.0_Dev/training/compendium_data/microarray_CV_ref_list.rds")
labels <- get_label_list(reference)

createSignaturesV1() <- function(reference, labels, quantile_interval = 0.0001, min_prob = 0.0001, max_genes = 100, min_genes = 5,
                               diff_vals, use_counts = FALSE, plot = TRUE){

  if (is.list(reference)) {
    for (ref in names(reference)) {
      ref_name <- "hpca" # Remove

      ref <- reference[[ref_name]]  # SE object
      ref_labels <- labels[[ref_name]]  # Cell types

      celltype_signatures_list <- list()


      for (ctoi in ref_labels) {
        overall_ctoi_up_genes <- c()
        overall_ctoi_down_genes <- c()

        for (celltype in ref_labels) {
          #ctoi <- "Astrocyte" # Remove
          #celltype <- "B cell" # Remove

          if (ctoi == celltype) {
            next
          }

          # Rank genes
          ctoi_ranked <- rankGenes(assays(ref[, ref$label.fine == ctoi])$logcounts)
          celltype_ranked <- rankGenes(assays(ref[, ref$label.fine == celltype])$logcounts)

          # Get the mean rank
          ctoi_ranked_avg <-  apply(ctoi_ranked, 1, mean)
          celltype_ranked_avg <- apply(celltype_ranked, 1, mean)

          # Get the difference in mean ranks
          rank_diff <- ctoi_ranked_avg - celltype_ranked_avg

          # Find up and down expressed genes in the CTOI
          up_genes <- c()
          down_genes <- c()
          while (length(up_genes) < min_genes && length(down_genes) < min_genes) {
            up_genes <- unique(c(up_genes, names(rank_diff[rank_diff > quantile(rank_diff, 1-min_prob)])))
            down_genes <- unique(c(down_genes, names(rank_diff[rank_diff < quantile(rank_diff, min_prob)])))
            min_prob <- min_prob + quantile_interval
          }

          # Keep up and down genes against all cell types
          overall_ctoi_up_genes <- c(overall_ctoi_up_genes, up_genes)
          overall_ctoi_down_genes <- c(overall_ctoi_down_genes, down_genes)
        }


        # min_prob <- min_prob - quantile_interval
        # Plot histogram of rank_diff quantiles
        # if (plot) {
        #   df <- data.frame(celltype = "T cell", rank_diff = rank_diff, genes = names(rank_diff))
        #   ggplot(df, aes(x=celltype, y=rank_diff, col=celltype)) +
        #     scale_y_continuous(breaks = c(quantile(rank_diff, min_prob), quantile(rank_diff, 1-min_prob))) +
        #     geom_point(position = position_jitter(seed = 1, width = 0.2), alpha = 0.2) +
        #     geom_violin() +
        #     geom_label_repel(data = df[df$genes %in% c(up_genes, down_genes),], aes(celltype, rank_diff, label = genes),
        #                     color = "black", size=3, box.padding = 0.5) +
        #     labs(y="Difference in gene rank", x="", col="") +
        #     theme_minimal(base_size = 11)
        # }


      }
    }
  }else{

    if (use_counts) {

      else{

      }
    }
  }
}


create







## Creating signatures with counts of genes (Florian) ----
createSignatures_flo <- function(expr, main_labels, probs, diff_vals) {
  celltype_signatures_hash <- hash()
  for(celltype in main_labels) {
    ## subset reference data
    in_ref <- assays(expr[, expr$label.main == celltype])$logcount
    out_ref <- assays(expr[, expr$label.main != celltype])$logcount
    ## get quantiles for each gene in both subsets:
    in_q <- apply(in_ref, 1, function(x) {quantile(x, probs, na.rm = TRUE)})
    out_q <- apply(out_ref, 1, function(x) {quantile(x, 1-probs, na.rm = TRUE)})
    gene_ids <- row.names(expr)
    #cf <- matrix(nrow = length(diff_vals)); len_probs <- length(probs)
    signature_list <- list()
    for(pidx in seq(length(probs))) {
      p <- paste0("prob_",probs[pidx])
      for(x in diff_vals) {
        d <- paste0("diff_",x)
        gene_vec <- gene_ids[as.vector(in_q[pidx, ] > (out_q[pidx, ] + x))]
        if(length(gene_vec) > 8 && length(gene_vec) < 200) {
          gs_name <- paste(celltype,d,p,sep = "+")
          gs <- GeneSet(gene_vec, setName = gs_name)
          signature_list[[gs_name]] <- gs
        }
      }
    }
    celltype_signatures_hash[[celltype]] <- signature_list
  }
  return(celltype_signatures_hash)
}



# Required libraries and functions:
library(singscore)
library(SummarizedExperiment)
library(GSEABase)

get_label_list <- function(ref_list) {
  lapply(ref_list, function(ref) {
    levels(as.factor(ref$label.fine))
  })
}




random_ref_list <- readRDS(file = "~/Documents/xCell2.0/xCell2.0_Dev/training/compendium_data/microarray_CV_ref_list.rds")
label_list <- get_label_list(random_ref_list)


## Creating signatures with ranks of genes (Florian) ----
create_signatures_16_01 <- function(random_ref_list, label_list, num_quantiles = 20, plotting = TRUE) {

  # Prob hold different cutoffs for quantiles
  probs = seq(from = 0.0001,to = 0.1,length.out = num_quantiles)

  csl <- lapply(names(random_ref_list), function(ref_name) {  # For each reference
    ref_name <- "hpca" # Remove

    ref <- random_ref_list[[ref_name]]  # SE object

    ref_labels <- label_list[[ref_name]]  # Cell types

    if(plotting) {par(mfrow = c(ceiling(length(ref_labels)/3), 3))}   # ?

    celltype_signatures_list <- list()

    for(celltype in ref_labels) {    # For each cell type
      celltype <- "T cell" # Remove

      ## ranks for genes in- and outside of samples for CTOI
      in_df <- rankGenes(assays(ref[, ref$label.fine == celltype])$logcounts)
      out_df <- rankGenes(assays(ref[, ref$label.fine != celltype])$logcounts)

      ## get average rank:
      in_df <- apply(in_df, 1, mean)
      out_df <- apply(out_df, 1, mean)

      ## get difference in ranks:
      rank_diff <- in_df - out_df

      ## get quantiles:
      q_diff <- quantile(rank_diff, c(probs, 1-probs[length(probs):1]))

      if(plotting) {hist(rank_diff, breaks = 50, main = paste0("Histogram of rank-difference: ", celltype)); abline(v = q_diff)}

      for(p in probs) {    # For each quantile
        p <- 0.000100000 # Remove

        ## get genes for corresponding quantile (up and downregulated):
        p_idx <- match(p,probs)

        genes_in_down <- names(rank_diff[rank_diff < q_diff[p_idx]])

        genes_in_up <- names(rank_diff[rank_diff > q_diff[length(q_diff)-p_idx+1]])

        if(length(genes_in_up) >= 5 && length(genes_in_up) < 100 && length(genes_in_down) >= 5 && length(genes_in_down) < 100) {
          gs_name <- paste(ref_name,celltype,p,sep = ",")
          gs_list <- list(
            "up"=GeneSet(genes_in_up[!is.na(genes_in_up)], setName = paste0(gs_name,",up")),
            "down"=GeneSet(genes_in_down[!is.na(genes_in_down)], setName = paste0(gs_name, ",down"))
          )
          celltype_signatures_list[[gs_name]] <- gs_list
        }
      }
    }


    celltype_signatures_list
  })


  if(plotting) {
    sig_sizes <- lapply(names(random_ref_list), function(ref_name) {
      ref <- random_ref_list[[ref_name]]
      ref_labels <- label_list[[ref_name]]
      sig_size_df <- matrix(ncol = length(ref_labels), nrow = length(probs)); colnames(sig_size_df) <- ref_labels; rownames(sig_size_df) <- probs
      for(celltype in ref_labels) {
        in_df <- rankGenes(assays(ref[, ref$label.fine == celltype])$logcounts)
        out_df <- rankGenes(assays(ref[, ref$label.fine != celltype])$logcounts)
        in_df <- apply(in_df, 1, mean); out_df <- apply(out_df, 1, mean)
        rank_diff <- in_df - out_df
        q_diff <- quantile(rank_diff, c(probs,1-probs[length(probs):1]))
        for(p in probs) {
          p_idx <- match(p,probs)
          genes_in_up <- names(rank_diff[rank_diff > q_diff[length(q_diff)-p_idx+1]])
          sig_size_df[p_idx,match(celltype, ref_labels)] <- length(genes_in_up)
        }
      }
      sig_size_df
    })
    names(sig_sizes) <- names(random_ref_list)
    for(ref_name in names(sig_sizes)) {
      print(pheatmap::pheatmap(sig_sizes[[ref_name]], cluster_cols = F, cluster_rows = F, angle_col = 45, display_numbers = TRUE,
                               main = paste("Number of genes for different values of p (reference:", ref_name,")")))
    }
  }
  names(csl) <- names(random_ref_list)
  return(csl)
}




create_signatures_out <- create_signatures_16_01(ref_list, lab_list, num_quantiles = 20, plotting = F)
