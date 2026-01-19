library(Seurat)
library(spacexr)
library(Matrix)

#' Run RCTD for a single Visium sample given a Seurat reference
#'
#' @param sobj                 Seurat object for the spatial sample
#' @param sc_ref               Seurat object for the single cell reference
#' @param assay                Assay name used in the spatial Seurat object
#' @param slot                 Slot name used in the spatial assay (for example "counts")
#' @param sc_ref_celltype_col  Column in sc_ref meta.data containing cell type labels
#' @param sc_ref_slot          Slot name used in the reference assay (for example "counts")
#' @param cores                Maximum number of cores for spacexr::create.RCTD
#'
#' @return Seurat object sobj with a new assay "rctd_full" and per cell type weights in meta.data
runRCTD_1_sample <- function(
  sobj,
  sc_ref,
  assay               = "SCT",
  slot                = "counts",
  sc_ref_celltype_col = "Celltype_Mal_sub",
  sc_ref_slot         = "counts",
  cores               = 2
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

  ## Create SpatialRNA object for RCTD query
  ## Uses the selected assay and slot from the spatial Seurat object
  counts <- Seurat::GetAssayData(sobj, assay = assay, slot = slot)

  ## Get spatial coordinates and restrict to x and y
  coords <- Seurat::GetTissueCoordinates(sobj)[, c("x", "y")]

  ## Total UMI counts per spot needed by spacexr
  nUMI <- Matrix::colSums(counts)

  ## Construct spacexr SpatialRNA object
  query <- spacexr::SpatialRNA(coords, counts, nUMI)

  ## Extract reference counts matrix from the reference Seurat object
  ref_counts <- Seurat::GetAssayData(sc_ref, slot = sc_ref_slot)

  ## Fetch cell type labels from metadata using a stable Seurat API
  ## FetchData returns a data frame with one column per requested variable
  df <- Seurat::FetchData(sc_ref, vars = sc_ref_celltype_col)
  cell_types <- df[[1]]
  names(cell_types) <- rownames(df)   # barcodes as names

  ## Coerce to factor as expected by spacexr
  cell_types <- as.factor(cell_types)

  ## Ensure that the cell type vector is aligned to the reference counts columns
  ## This guarantees that each column in ref_counts has a label
  if (!all(colnames(ref_counts) %in% names(cell_types))) {
    stop(
      "Not all reference barcodes have a cell type label. ",
      "Check that sc_ref_celltype_col is complete and matches column names of ref_counts."
    )
  }

  ## Reorder and subset cell_types to match ref_counts columns exactly
  cell_types <- cell_types[colnames(ref_counts)]

  ## Total UMI per reference cell
  nUMI_ref <- Matrix::colSums(ref_counts)

  ## Construct spacexr Reference object
  reference <- spacexr::Reference(ref_counts, cell_types, nUMI_ref)

  ## Run RCTD using spacexr
  ## create.RCTD builds the model object, run.RCTD performs the deconvolution
  rctd      <- spacexr::create.RCTD(query, reference, max_cores = cores)
  rctd_full <- spacexr::run.RCTD(rctd, doublet_mode = "full")

  ## Extract and normalize RCTD weights
  ## weights is a cell type by spot matrix, normalize_weights converts to probabilities per spot
  weights      <- rctd_full@results$weights
  norm_weights <- spacexr::normalize_weights(weights)

  ## Add RCTD results as a new assay on the spatial Seurat object
  ## Note the transpose so that rows become features (cell types) and columns spots
  sobj[["rctd_full"]] <- Seurat::CreateAssayObject(
    data = t(as.matrix(norm_weights))
  )

  ## Ensure the assay key is set to a consistent prefix for downstream usage
  if (is.null(sobj@assays$rctd_full@key) || length(sobj@assays$rctd_full@key) == 0) {
    sobj@assays$rctd_full@key <- "rctd_full_"
  }

  ## Add one metadata column per cell type to the spatial Seurat object
  ## Column names will match the row names of norm_weights (cell types)
  sobj <- Seurat::AddMetaData(
    object   = sobj,
    metadata = norm_weights
  )

  ## Return the enriched spatial Seurat object
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
