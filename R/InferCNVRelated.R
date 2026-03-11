library(Seurat)
library(dplyr)

#' Add inferCNV HMM-based CNV scores to a Seurat V5 object
#'
#' Reads inferCNV HMM output files, maps predicted CNV states from cell groups
#' back to individual cells, computes per-cell CNV summary scores, and adds
#' those scores to the Seurat object's meta.data.
#'
#' @param sobj                Seurat object to annotate
#' @param cell_groupings_path Path to inferCNV cell groupings file. Must contain
#'                            columns "cell_group_name" and "cell"
#' @param pred_cnv_genes_path Path to inferCNV predicted CNV genes file. Must
#'                            contain columns "cell_group_name" and "state"
#' @param neutral_state       HMM state treated as copy-number neutral
#' @param amp_states          Vector of HMM states treated as amplifications
#' @param del_states          Vector of HMM states treated as deletions
#' @param prefix              Optional prefix added to output metadata column names
#'
#' @return Seurat object sobj with added meta.data columns:
#'         "total_cnv", "cnv_amp_score", and "cnv_del_score"
add_infercnv_hmm_scores_to_seurat <- function(
  sobj,
  cell_groupings_path,
  pred_cnv_genes_path,
  neutral_state = 3,
  amp_states = c(4, 5, 6),
  del_states = c(1, 2),
  prefix = ""
) {
  stopifnot(file.exists(cell_groupings_path), file.exists(pred_cnv_genes_path))

  suppressPackageStartupMessages({
    require(dplyr)
    require(Seurat)
  })

  cell_groups <- read.table(cell_groupings_path, header = TRUE, stringsAsFactors = FALSE)
  cnv_genes   <- read.table(pred_cnv_genes_path, header = TRUE, stringsAsFactors = FALSE)

  required_cg <- c("cell_group_name", "cell")
  required_cnv <- c("cell_group_name", "state")

  missing_cg <- setdiff(required_cg, colnames(cell_groups))
  missing_cnv <- setdiff(required_cnv, colnames(cnv_genes))

  if (length(missing_cg) > 0) stop("cell_groupings is missing columns: ", paste(missing_cg, collapse = ", "))
  if (length(missing_cnv) > 0) stop("pred_cnv_genes is missing columns: ", paste(missing_cnv, collapse = ", "))

  cnv_data <- cnv_genes %>%
    left_join(cell_groups, by = "cell_group_name")

  if (!"cell" %in% colnames(cnv_data)) {
    stop("After join, column 'cell' is missing. Check that cell_groupings has 'cell' and 'cell_group_name'.")
  }

  cnv_scores <- cnv_data %>%
    group_by(cell) %>%
    summarise(
      total_cnv      = sum(state != neutral_state, na.rm = TRUE),
      cnv_amp_score  = sum(state %in% amp_states,  na.rm = TRUE),
      cnv_del_score  = sum(state %in% del_states,  na.rm = TRUE),
      .groups = "drop"
    )

  add_named <- function(vec, cells, colname) {
    names(vec) <- cells
    Seurat::AddMetaData(sobj, metadata = vec, col.name = colname)
  }

  sobj <- add_named(cnv_scores$total_cnv,     cnv_scores$cell, paste0(prefix, "total_cnv"))
  sobj <- add_named(cnv_scores$cnv_amp_score, cnv_scores$cell, paste0(prefix, "cnv_amp_score"))
  sobj <- add_named(cnv_scores$cnv_del_score, cnv_scores$cell, paste0(prefix, "cnv_del_score"))

  overlap <- sum(cnv_scores$cell %in% colnames(sobj))
  message("Added CNV metadata. Cells matched to Seurat object: ", overlap, " / ", nrow(cnv_scores))

  return(sobj)
}
