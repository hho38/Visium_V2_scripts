#' Running SpotClean on Multiple Samples 
#' Original Paper : https://www.nature.com/articles/s41467-022-30587-y

#' Load Visium data, run SpotClean, optionally visualize, and return a Seurat object
#'
#' @description
#' Convenience wrapper that
#' 1) validates a Space Ranger output directory,
#' 2) loads the raw counts and spatial metadata,
#' 3) builds a slide object,
#' 4) runs decontamination with `spotclean`,
#' 5) optionally visualizes the estimated contamination rate as a heatmap,
#' 6) converts the result to a Seurat object,
#' 7) returns all key objects in a single list.
#'
#' @param data_dir Character. Path to a Space Ranger run directory that contains
#'   `raw_feature_bc_matrix.h5` and a `Spatial` subfolder with
#'   `tissue_positions.csv`, `tissue_lowres_image.png`, and `scalefactors_json.json`.
#'   Default is `"./spaceranger_data"`.
#'
#' @param visualize Logical. If `TRUE`, draw a heatmap of `contamination_rate`
#'   stored in `decont_obj@metadata`. Default is `TRUE`.
#'
#' @param legend_title Character. Title for the heatmap legend. Default is `"contamination rate"`.
#'
#' @param legend_range Numeric length two. Limits for the heatmap legend range.
#'   Default is `c(0, 1)`.
#'
#' @param spotclean_args List. Named arguments forwarded to `spotclean` via `do.call`.
#'   Use this to tune the decontamination step.
#'
#' @param seurat_args List. Named arguments forwarded to `convertToSeurat` via `do.call`.
#'   Use this to control Seurat conversion and image attachment.
#'
#' @details
#' The function checks that the following functions exist in the current session:
#' `read10xRawH5`, `read10xSlide`, `createSlide`, `spotclean`, `visualizeHeatmap`,
#' and `convertToSeurat`. Load the packages that provide these before calling.
#'
#' Required files are resolved from `data_dir`:
#' - `raw_feature_bc_matrix.h5`
#' - `Spatial/tissue_positions.csv`
#' - `Spatial/tissue_lowres_image.png`
#' - `Spatial/scalefactors_json.json`
#'
#' @return
#' A list with three elements:
#' - `slide_obj`: the slide object created from raw counts and spatial metadata
#' - `decont_obj`: the decontaminated slide object returned by `spotclean`
#' - `seurat`: a Seurat object created by `convertToSeurat`
#'
#' @examples
#' \dontrun{
#' # Minimal run with defaults
#' out <- load_and_SpotClean(data_dir = "path/to/spaceranger")
#'
#' # Customize SpotClean and Seurat steps
#' out <- load_and_SpotClean(
#'   data_dir = "path/to/spaceranger",
#'   visualize = TRUE,
#'   legend_title = "estimated contamination",
#'   legend_range = c(0, 0.5),
#'   spotclean_args = list(max_iters = 200, seed = 123),
#'   seurat_args = list(project = "Sample1", assay = "Spatial")
#' )
#'
#' # Access returned objects
#' out$seurat
#' out$decont_obj@metadata$contamination_rate
#' }
#'
#' @seealso
#' `spotclean`, `visualizeHeatmap`, `convertToSeurat`
#'
#' @export
#' @importFrom spotclean read10xRawH5 read10xSlide createSlide spotclean visualizeHeatmap convertToSeurat
#' @import Seurat


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