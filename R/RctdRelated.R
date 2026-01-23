library(Seurat)
library(spacexr)
library(Matrix)

#' Run RCTD for a single Visium sample given a Seurat reference
#'
#' @param sobj                 Seurat object for the spatial sample
#' @param sc_ref               Seurat object for the single cell reference
#' @param assay                Assay name used in the spatial Seurat object (Spatial, SCT, Ect)
#' @param slot                 Slot name used in the spatial assay (for example "counts" or "data")
#' @param sc_ref_celltype_col  Column in sc_ref meta.data containing cell type labels
#' @param sc_ref_slot          Slot name used in the reference assay (for example "counts")
#' @param cores                Maximum number of cores for spacexr::create.RCTD
#' @param min_cells_per_type   Minimum number of cells per celltype in reference. If lower, celltype won't be used. 
#'
#' @return Seurat object sobj with a new assay "rctd_full" and per cell type weights in meta.data
runRCTD_1_sample <- function(
    sobj,
    sc_ref,
    assay               = "SCT",
    slot                = "counts",
    sc_ref_celltype_col = "Celltype_col",
    sc_ref_slot         = "counts",
    cores               = 2,
    min_cells_per_type  = 25
) {
  ## Basic input validation to fail early with clear messages
  if (!inherits(sobj, "Seurat")) {
    stop("Argument 'sobj' must be a Seurat object.")
  }
  if (!inherits(sc_ref, "Seurat")) {
    stop("Argument 'sc_ref' must be a Seurat object.")
  }
  if (!sc_ref_celltype_col %in% colnames(sc_ref@meta.data)) {
    stop(
      "Column '", sc_ref_celltype_col,
      "' not found in sc_ref@meta.data. Check sc_ref_celltype_col."
    )
  }
  if (!is.numeric(min_cells_per_type) || length(min_cells_per_type) != 1 || min_cells_per_type < 1) {
    stop("min_cells_per_type must be a single numeric value >= 1.")
  }
  
  ## Create SpatialRNA object for RCTD query
  counts <- Seurat::GetAssayData(sobj, assay = assay, slot = slot)
  coords <- Seurat::GetTissueCoordinates(sobj)[, c("x", "y")]
  nUMI   <- Matrix::colSums(counts)
  query  <- spacexr::SpatialRNA(coords, counts, nUMI)
  
  ## Extract reference counts matrix from the reference Seurat object
  ref_counts <- Seurat::GetAssayData(sc_ref, slot = sc_ref_slot)
  
  ## Fetch cell type labels
  df <- Seurat::FetchData(sc_ref, vars = sc_ref_celltype_col)
  cell_types <- df[[1]]
  names(cell_types) <- rownames(df)
  cell_types <- as.factor(cell_types)
  
  ## Ensure labels cover reference barcodes
  if (!all(colnames(ref_counts) %in% names(cell_types))) {
    stop(
      "Not all reference barcodes have a cell type label. ",
      "Check that sc_ref_celltype_col is complete and matches column names of ref_counts."
    )
  }
  cell_types <- cell_types[colnames(ref_counts)]
  
  ## Filter out rare cell types below threshold
  ct_counts <- table(cell_types)
  keep_types <- names(ct_counts)[ct_counts >= min_cells_per_type]
  if (length(keep_types) < 1) {
    stop(
      "After filtering, no cell types remain with >= ", min_cells_per_type,
      " cells. Lower min_cells_per_type or check sc_ref_celltype_col."
    )
  }
  
  keep_cells <- cell_types %in% keep_types
  if (any(!keep_cells)) {
    message(
      "Filtering reference cell types: kept ", sum(keep_cells), " of ", length(keep_cells),
      " cells across ", length(keep_types), " cell types (min ", min_cells_per_type, " cells per type)."
    )
  }
  
  ref_counts <- ref_counts[, keep_cells, drop = FALSE]
  cell_types <- droplevels(cell_types[keep_cells])
  
  ## Total UMI per reference cell
  nUMI_ref <- Matrix::colSums(ref_counts)
  
  ## Construct spacexr Reference object
  reference <- spacexr::Reference(ref_counts, cell_types, nUMI_ref)
  
  ## Run RCTD
  rctd      <- spacexr::create.RCTD(query, reference, max_cores = cores)
  rctd_full <- spacexr::run.RCTD(rctd, doublet_mode = "full")
  
  ## Extract and normalize weights
  weights      <- rctd_full@results$weights
  norm_weights <- spacexr::normalize_weights(weights)
  
  ## Add as assay
  sobj[["rctd_full"]] <- Seurat::CreateAssayObject(
    data = t(as.matrix(norm_weights))
  )
  
  if (is.null(sobj@assays$rctd_full@key) || length(sobj@assays$rctd_full@key) == 0) {
    sobj@assays$rctd_full@key <- "rctd_full_"
  }
  
  ## Add per cell type weights to metadata
  sobj <- Seurat::AddMetaData(object = sobj, metadata = norm_weights)
  
  return(sobj)
}


# Pick, for every spot, the cell-type column that has the highest proportion.
# ─────────────────────────────────────────────────────────────────────────────
# object        : a Seurat object
# labels_col    : name of the new metadata column that will hold the chosen label
# celltype_cols : character vector of metadata columns that contain the
#                 per-cell-type prediction probabilities (one column per type)

assignLabels <- function(object,
                         labels_col    = "majority_celltype",
                         celltype_cols = NULL) {

  ## sanity checks -------------------------------------------------------------
  if (is.null(celltype_cols) || length(celltype_cols) == 0)
    stop("'celltype_cols' must be a non-empty character vector")

  missing <- setdiff(celltype_cols, colnames(object@meta.data))
  if (length(missing))
    stop("These columns are not in the metadata: ",
         paste(missing, collapse = ", "))

  ## pull the probability block -----------------------------------------------
  preds <- object@meta.data[, celltype_cols, drop = FALSE]

  if (!all(vapply(preds, is.numeric, logical(1))))
    stop("All columns in 'celltype_cols' must be numeric probabilities")

  ## choose the label with the highest probability for each spot --------------
  max_idx <- apply(preds, 1, which.max)             # index of the winner
  labels  <- celltype_cols[max_idx]                 # corresponding names

  ## write back to the object --------------------------------------------------
  object[[labels_col]] <- factor(labels, levels = celltype_cols)
  Idents(object)      <- labels_col

  return(object)
}
