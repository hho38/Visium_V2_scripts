#' Running SpotClean on Multiple Samples 


load_and_SpotClean <- function(
  data_dir = "./spaceranger_data",
  visualize = TRUE,
  legend_title = "contamination rate",
  legend_range = c(0, 1),
  spotclean_args = list(),
  seurat_args = list()
) {
  # Build paths
  spatial_dir       <- file.path(data_dir, "raw_feature_bc_matrix.h5")
  spatial_slide_dir <- file.path(data_dir, "Spatial", "tissue_positions.csv")
  spatial_img_dir   <- file.path(data_dir, "Spatial", "tissue_lowres_image.png")
  scale_factor_dir  <- file.path(data_dir, "Spatial", "scalefactors_json.json")

  # Basic file checks
  paths <- c(
    spatial_dir       = spatial_dir,
    spatial_slide_dir = spatial_slide_dir,
    spatial_img_dir   = spatial_img_dir,
    scale_factor_dir  = scale_factor_dir
  )
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "File not found for: ",
      paste(missing, collapse = ", "),
      " inside data_dir = ", normalizePath(data_dir)
    )
  }

  # Check that required functions exist in the session
  needed_funs <- c(
    "read10xRawH5", "read10xSlide", "createSlide",
    "spotclean", "visualizeHeatmap", "convertToSeurat"
  )
  missing_funs <- needed_funs[!vapply(needed_funs, function(f) exists(f, mode = "function"), logical(1))]
  if (length(missing_funs) > 0) {
    stop(
      "Missing functions in the current session: ",
      paste(missing_funs, collapse = ", "),
      ". Load the packages that provide them before calling this function."
    )
  }

  # Load data
  bc_raw   <- read10xRawH5(spatial_dir)
  bc_slide <- read10xSlide(spatial_slide_dir, spatial_img_dir, scale_factor_dir)
  slide_obj <- createSlide(bc_raw, bc_slide)

  # Decontaminate
  decont_obj <- do.call(
    "spotclean",
    c(list(slide_obj = slide_obj), spotclean_args)
  )

  # Optional visualization of contamination rate
  if (isTRUE(visualize)) {
    visualizeHeatmap(
      decont_obj,
      decont_obj@metadata$contamination_rate,
      logged = FALSE,
      legend_title = legend_title,
      legend_range = legend_range
    )
  }

  # Convert to Seurat
  sobj <- do.call(
    "convertToSeurat",
    c(list(decont_obj, image_dir = file.path(data_dir, "Spatial")), seurat_args)
  )

  # Return a handy bundle
  return(list(
    slide_obj  = slide_obj,
    decont_obj = decont_obj,
    seurat     = sobj
  ))
}