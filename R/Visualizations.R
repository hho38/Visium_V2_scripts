#' Plot Hallmark spatial signature scores with a shared color scale across samples as default Seurat will skew charts if pathway name is too long 
#'
#' This function generates spatial feature plots for one or more Hallmark signature
#' score features stored in a Seurat object. For each feature, the function:
#' \itemize{
#'   \item Computes global display limits across the full Seurat object using the
#'   1st and 99th percentiles.
#'   \item Creates one spatial plot per image/sample with identical color scaling.
#'   \item Combines the plots into a single figure with one shared legend.
#' }
#'
#' This is useful when comparing the same pathway or signature across multiple
#' Visium samples, because all panels use the same value range and palette.
#'
#' @param sobj A Seurat object containing spatial data and the requested features.
#' @param features A character vector of feature names to plot. Each feature must
#'   be retrievable with [Seurat::FetchData()].
#' @param ncol Number of columns used when arranging the per-sample plots.
#'   Default is `4`.
#' @param pt.size.factor Numeric scaling factor for spot size in
#'   [Seurat::SpatialFeaturePlot()]. Default is `7`.
#'
#' @return Invisibly returns `NULL`. The function is called for its side effect of
#'   printing one combined plot per feature.
#'
#' @details
#' The plotting range is defined independently for each feature using the 1st and
#' 99th percentiles of that feature across the entire object. Values outside this
#' range are squished into the displayed limits using [scales::squish()].
#'
#' The legend title is cleaned by removing a leading `"HALLMARK_"` prefix and a
#' trailing `"_UCell"` suffix.
#'
#' @examples
#' \dontrun{
#' plot_hallmark_spatial_signatures_shared_scale(
#'   sobj = all_combined_sobj,
#'   features = c("HALLMARK_HYPOXIA_UCell", "HALLMARK_IL6_JAK_STAT3_SIGNALING_UCell"),
#'   ncol = 4,
#'   pt.size.factor = 7
#' )
#' }
#'
#' @importFrom Seurat FetchData SpatialFeaturePlot
#' @importFrom ggplot2 scale_fill_gradientn theme_minimal theme element_blank
#' @importFrom ggplot2 element_text margin
#' @importFrom patchwork wrap_plots
#' @importFrom grid unit
#' @importFrom scales squish
#' @export
plot_hallmark_spatial_signatures_shared_scale <- function(
  sobj,
  features,
  ncol = 4,
  pt.size.factor = 7
) {
  # ---- Input validation ------------------------------------------------------

  if (!inherits(sobj, "Seurat")) {
    stop("`sobj` must be a Seurat object.", call. = FALSE)
  }

  if (!is.character(features) || length(features) == 0) {
    stop("`features` must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) || ncol < 1) {
    stop("`ncol` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(pt.size.factor) || length(pt.size.factor) != 1 || is.na(pt.size.factor) || pt.size.factor <= 0) {
    stop("`pt.size.factor` must be a single positive number.", call. = FALSE)
  }

  # ---- Iterate through requested features -----------------------------------

  for (f in features) {
    # Confirm that the feature exists and can be fetched from the object.
    vals <- tryCatch(
      FetchData(sobj, vars = f, layer = "data")[[1]],
      error = function(e) {
        stop(sprintf("Feature '%s' could not be fetched from `sobj`.", f), call. = FALSE)
      }
    )

    # Ensure the feature contains usable numeric values.
    if (!is.numeric(vals)) {
      stop(sprintf("Feature '%s' is not numeric and cannot be plotted on a continuous scale.", f), call. = FALSE)
    }

    if (all(is.na(vals))) {
      warning(sprintf("Feature '%s' contains only NA values. Skipping.", f), call. = FALSE)
      next
    }

    # Use the 1st and 99th percentiles to reduce the influence of extreme outliers
    # while preserving a shared scale across all samples/images in the object.
    lims <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)

    # Build a cleaner legend title by removing common Hallmark/UCell wrappers.
    clean_name <- gsub("_UCell$", "", gsub("^HALLMARK_", "", f))
    legend_title <- clean_name

    # Generate one SpatialFeaturePlot per image/sample without combining them yet.
    # Using `combine = FALSE` allows us to manually apply identical styling and
    # collect a single shared legend afterward.
    plist <- SpatialFeaturePlot(
      sobj,
      features = f,
      combine = FALSE,
      pt.size.factor = pt.size.factor,
      min.cutoff = lims[1],
      max.cutoff = lims[2]
    )

    # Apply the exact same continuous color scale to every plot so that values
    # are directly comparable across samples.
    plist <- lapply(plist, function(pp) {
      pp +
        scale_fill_gradientn(
          colours = c(
            "#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
            "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
          ),
          limits = lims,
          oob = squish,
          name = legend_title
        ) +
        theme_minimal(base_size = 9) +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.7, "cm"),
          plot.margin = margin(4, 4, 4, 4)
        )
    })

    # Combine per-sample plots into a single figure and collect all legends into
    # one shared legend positioned on the right.
    p <- wrap_plots(plist, ncol = ncol, guides = "collect") &
      theme(legend.position = "right")

    # Print the combined figure for the current feature.
    print(p)
  }

  invisible(NULL)
}
