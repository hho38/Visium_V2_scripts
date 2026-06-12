#' Spatial blended feature plot for two features
#'
#' Creates spatial plots for two features and a blended overlap plot. The first
#' feature is mapped along the x-axis of the blend color grid, and the second
#' feature is mapped along the y-axis. The corner colors define the appearance
#' of low/low, high feature 1 only, high feature 2 only, and high/high overlap.
#'
#' This function supports percentile clipping, optional binned colors, and
#' optional shared color scaling across multiple spatial images/samples.
#'
#' @param object A Seurat object containing spatial image data.
#'
#' @param features A character vector of length 2 giving the two features to
#'   blend. These can be genes or metadata columns. `features[1]` is treated as
#'   the x feature and maps from `bottom_left` to `bottom_right`.
#'   `features[2]` is treated as the y feature and maps from `bottom_left` to
#'   `top_left`.
#'
#' @param combine Logical. If `TRUE`, returns one combined patchwork plot. If
#'   `FALSE`, returns a nested list of plots for each image/sample. Default is
#'   `TRUE`.
#'
#' @param feature_1_alt_name Optional character string used as the x-axis label
#'   in the blend legend for `features[1]`. If `NULL`, the original feature name
#'   is used.
#'
#' @param feature_2_alt_name Optional character string used as the y-axis label
#'   in the blend legend for `features[2]`. If `NULL`, the original feature name
#'   is used.
#'
#' @param assay Optional assay name to set as the default assay before plotting.
#'   For example, `"SCT"` or `"Spatial"`. If `NULL`, the current default assay is
#'   used.
#'
#' @param bottom_left Color for spots with low values for both features. This is
#'   also the low end of each single-feature color scale. Default is black
#'   (`"#000000"`).
#'
#' @param bottom_right Color for spots with high values for feature 1 and low
#'   values for feature 2. This is the high color for `features[1]`. Default is
#'   red (`"#FF0000"`).
#'
#' @param top_left Color for spots with low values for feature 1 and high values
#'   for feature 2. This is the high color for `features[2]`. Default is green
#'   (`"#00FF00"`).
#'
#' @param top_right Color for spots with high values for both features. This is
#'   the strongest overlap color. Default is yellow (`"#FFFF00"`).
#'
#' @param x_bottom_percentile Numeric value from 0 to 100. Lower percentile
#'   cutoff for feature 1. Values below this percentile are clipped to the cutoff
#'   before color mapping. Default is `0`.
#'
#' @param x_top_percentile Numeric value from 0 to 100. Upper percentile cutoff
#'   for feature 1. Values above this percentile are clipped to the cutoff before
#'   color mapping. For example, `90` treats all values above the 90th percentile
#'   as the maximum visible feature-1 value. Default is `100`.
#'
#' @param y_bottom_percentile Numeric value from 0 to 100. Lower percentile
#'   cutoff for feature 2. Values below this percentile are clipped to the cutoff
#'   before color mapping. Default is `0`.
#'
#' @param y_top_percentile Numeric value from 0 to 100. Upper percentile cutoff
#'   for feature 2. Values above this percentile are clipped to the cutoff before
#'   color mapping. For example, `90` treats all values above the 90th percentile
#'   as the maximum visible feature-2 value. Default is `100`.
#'
#' @param uniform_color_scale Logical. If `TRUE`, percentile limits and bin
#'   breaks are calculated once using all spatial spots across all images, so the
#'   same numeric value maps to the same color or bin in every sample. If
#'   `FALSE`, each image/sample is scaled independently. Default is `TRUE`.
#'
#' @param n_bins Integer. Controls whether colors are continuous or binned.
#'   Use `n_bins = 0` for continuous colors. Use `n_bins > 1` for discrete
#'   binned colors. For example, `n_bins = 4` creates 4 bins per feature and
#'   16 possible blended colors. `n_bins = 1` is not meaningful and will return
#'   an error. Default is `0`.
#'
#' @param use_seurat_backend Logical. If `TRUE`, uses Seurat's `FeaturePlot`
#'   blend backend after creating a spatial dimensional reduction. If `FALSE`,
#'   uses the custom color-grid backend. The custom backend is recommended when
#'   using `top_right`, percentile clipping, binned colors, or
#'   `uniform_color_scale`. Default is `FALSE`.
#'
#' @param fp_extra_arguments Named list of extra arguments passed to Seurat's
#'   `FeaturePlot()` when `use_seurat_backend = TRUE`. If this list is not empty,
#'   `use_seurat_backend` is automatically set to `TRUE`.
#'
#' @param sfp_extra_arguments Named list of extra arguments passed to
#'   `SpatialDimPlot()`. Common examples include `pt.size.factor`,
#'   `image.alpha`, and other spatial plotting options.
#'
#' @return If `combine = TRUE`, returns a patchwork object containing the two
#'   single-feature spatial plots, the blended spatial plot, and the blend
#'   legend for each image. If `combine = FALSE`, returns a nested list of plots
#'   by image/sample.
#'
#' @details
#' The blend color grid is interpreted as:
#'
#' \itemize{
#'   \item `bottom_left`: low feature 1, low feature 2
#'   \item `bottom_right`: high feature 1, low feature 2
#'   \item `top_left`: low feature 1, high feature 2
#'   \item `top_right`: high feature 1, high feature 2
#' }
#'
#' Percentile cutoffs do not remove spots. Instead, values below the bottom
#' percentile are clipped to the bottom cutoff, and values above the top
#' percentile are clipped to the top cutoff.
#'
#' When `uniform_color_scale = TRUE`, these cutoffs are computed globally across
#' all spatial images in the object. This means a value such as `0.5` is mapped
#' to the same color or bin regardless of which sample/image it comes from.
#'
#' @examples
#' # Continuous blended colors with shared scaling across all samples
#' SpatialFeaturePlotBlend(
#'   object = all_combined_sobj,
#'   features = c("Macrophage", "Malignant"),
#'   top_left = "blue",
#'   bottom_right = "red",
#'   bottom_left = "white",
#'   top_right = "black",
#'   assay = "SCT",
#'   x_bottom_percentile = 1,
#'   x_top_percentile = 90,
#'   y_bottom_percentile = 1,
#'   y_top_percentile = 90,
#'   uniform_color_scale = TRUE,
#'   n_bins = 0,
#'   sfp_extra_arguments = list(pt.size.factor = 9)
#' )
#'
#' # Binned blended colors with shared bin cutoffs across all samples
#' SpatialFeaturePlotBlend(
#'   object = all_combined_sobj,
#'   features = c("Macrophage", "Malignant"),
#'   top_left = "blue",
#'   bottom_right = "red",
#'   bottom_left = "white",
#'   top_right = "black",
#'   assay = "SCT",
#'   x_bottom_percentile = 1,
#'   x_top_percentile = 90,
#'   y_bottom_percentile = 1,
#'   y_top_percentile = 90,
#'   uniform_color_scale = TRUE,
#'   n_bins = 4,
#'   sfp_extra_arguments = list(pt.size.factor = 9)
#' )
#'
#' @import Seurat
#' @import ggplot2
#' @import patchwork
#'
#' @export
SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL,
                                    assay = NULL,
                                    bottom_left = "lightgray",
                                    bottom_right = "red",
                                    top_left = "blue",
                                    top_right = "black",
                                    x_bottom_percentile = 0,
                                    x_top_percentile = 100,
                                    y_bottom_percentile = 0,
                                    y_top_percentile = 100,
                                    uniform_color_scale = TRUE,
                                    n_bins = 0,
                                    use_seurat_backend = FALSE,
                                    fp_extra_arguments = list(),
                                    sfp_extra_arguments = list()) {

    # Generate a 2D color grid.
    # x-axis: feature 1, from bottom_left to bottom_right
    # y-axis: feature 2, from bottom_left to top_left
    # high feature 1 + high feature 2 becomes top_right.
    gen_color_grid <- function(side_length, bottom_left, bottom_right,
                               top_left, top_right) {

        grad_gen <- function(start, end, n = side_length) {
            colorRampPalette(c(start, end))(n)
        }

        bottom_grad <- grad_gen(bottom_left, bottom_right)
        top_grad <- grad_gen(top_left, top_right)

        col_grid <- matrix(
            NA_character_,
            nrow = side_length,
            ncol = side_length
        )

        for (x_i in seq_len(side_length)) {
            col_grid[x_i, ] <- grad_gen(
                start = bottom_grad[x_i],
                end = top_grad[x_i]
            )
        }

        return(col_grid)
    }

    validate_percentiles <- function(bottom_percentile, top_percentile,
                                     feature_name = "feature") {

        if (bottom_percentile < 0 || bottom_percentile > 100 ||
            top_percentile < 0 || top_percentile > 100) {
            stop("Percentiles must be between 0 and 100.")
        }

        if (bottom_percentile >= top_percentile) {
            stop(paste0(
                "Bottom percentile must be lower than top percentile for ",
                feature_name,
                "."
            ))
        }
    }

    get_percentile_limits <- function(x, bottom_percentile, top_percentile,
                                      feature_name = "feature") {

        validate_percentiles(
            bottom_percentile = bottom_percentile,
            top_percentile = top_percentile,
            feature_name = feature_name
        )

        qs <- quantile(
            x,
            probs = c(bottom_percentile, top_percentile) / 100,
            na.rm = TRUE
        )

        if (any(is.na(qs))) {
            stop(paste0(
                "Could not compute percentile limits for ",
                feature_name,
                ". Check for missing or non-numeric values."
            ))
        }

        return(qs)
    }

    make_bin_breaks <- function(limits, n_bins) {
        seq(
            from = limits[1],
            to = limits[2],
            length.out = n_bins + 1
        )
    }

    # Map feature values to color-grid indices.
    # n_bins = 0 gives continuous scaling.
    # n_bins > 1 gives discrete bins.
    map_to_color_grid <- function(x, limits, side_length, n_bins = 0,
                                  bin_breaks = NULL) {

        x_clipped <- pmin(pmax(x, limits[1]), limits[2])

        if (limits[2] == limits[1]) {
            return(rep(1, length(x)))
        }

        if (n_bins > 1) {

            if (is.null(bin_breaks)) {
                bin_breaks <- make_bin_breaks(
                    limits = limits,
                    n_bins = n_bins
                )
            }

            grid_index <- findInterval(
                x = x_clipped,
                vec = bin_breaks,
                rightmost.closed = TRUE,
                all.inside = TRUE
            )

            grid_index <- pmin(pmax(grid_index, 1), n_bins)

        } else {

            x_scaled <- (x_clipped - limits[1]) / (limits[2] - limits[1])

            grid_index <- round((side_length - 1) * x_scaled) + 1
            grid_index <- pmin(pmax(grid_index, 1), side_length)
        }

        grid_index[is.na(grid_index)] <- 1

        return(grid_index)
    }

    custom_color_SpatialDimPlot <- function(cells_obj, image_name,
                                            new_md_column_name,
                                            colors_per_spot, ...) {

        cells_obj[[new_md_column_name]] <- colors_per_spot
        names(colors_per_spot) <- as.character(colors_per_spot)

        p <- SpatialDimPlot(
            object = cells_obj,
            group.by = new_md_column_name,
            cols = colors_per_spot,
            images = image_name,
            ...
        ) +
            ggtitle(new_md_column_name) +
            blend_plot_theme

        return(p)
    }

    extract_colors_from_ggplot <- function(p) {
        built <- ggplot_build(p)$data[[1]]

        if (!is.na(built[1, "fill"])) {
            col_to_use <- "fill"
        } else {
            col_to_use <- "colour"
        }

        return(built[, col_to_use])
    }

    if (length(features) != 2) {
        stop(paste0(
            "Incorrect number of features. Requires two features, received ",
            length(features),
            "."
        ))
    }

    if (n_bins < 0) {
        stop("n_bins must be 0 or greater. Use n_bins = 0 for continuous colors.")
    }

    if (n_bins == 1) {
        stop(
            paste0(
                "n_bins = 1 is not meaningful. ",
                "Use n_bins = 0 for continuous colors, or n_bins > 1 for binned colors."
            )
        )
    }

    use_binned_colors <- n_bins > 1

    if (!is.null(assay)) {
        DefaultAssay(object) <- assay
    }

    if (length(fp_extra_arguments) > 0) {
        use_seurat_backend <- TRUE
    }

    if (uniform_color_scale && use_seurat_backend) {
        warning(paste(
            "uniform_color_scale is only applied in the custom backend.",
            "Set use_seurat_backend = FALSE to use shared scaling."
        ))
    }

    blend_plot_theme <- theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
    )

    # -------------------------------------------------------------------------
    # Compute global limits and global bin breaks once if uniform scaling is on.
    # -------------------------------------------------------------------------

    global_x_lims <- NULL
    global_y_lims <- NULL
    global_x_breaks <- NULL
    global_y_breaks <- NULL

    if (!use_seurat_backend && uniform_color_scale) {

        image_cells <- unique(unlist(lapply(Images(object), function(img) {
            Seurat:::CellsByImage(
                object = object,
                images = img,
                unlist = TRUE
            )
        })))

        global_dat <- FetchData(
            object = object[, image_cells],
            vars = features
        )

        global_x_lims <- get_percentile_limits(
            x = global_dat[[features[1]]],
            bottom_percentile = x_bottom_percentile,
            top_percentile = x_top_percentile,
            feature_name = features[1]
        )

        global_y_lims <- get_percentile_limits(
            x = global_dat[[features[2]]],
            bottom_percentile = y_bottom_percentile,
            top_percentile = y_top_percentile,
            feature_name = features[2]
        )

        if (use_binned_colors) {
            global_x_breaks <- make_bin_breaks(
                limits = global_x_lims,
                n_bins = n_bins
            )

            global_y_breaks <- make_bin_breaks(
                limits = global_y_lims,
                n_bins = n_bins
            )
        }
    }

    plot_list_outer <- list()

    for (i in Images(object)) {

        cell_barcodes <- Seurat:::CellsByImage(
            object = object,
            images = i,
            unlist = TRUE
        )

        cells_obj_sub <- object[, cell_barcodes]

        images_sub_list <- list(object[[i]])
        names(images_sub_list) <- i
        cells_obj_sub@images <- images_sub_list

        if (!use_seurat_backend) {

            dat <- FetchData(
                object = cells_obj_sub,
                vars = features
            )

            side_length <- ifelse(use_binned_colors, n_bins, 100)

            col_grid <- gen_color_grid(
                side_length = side_length,
                bottom_left = bottom_left,
                bottom_right = bottom_right,
                top_left = top_left,
                top_right = top_right
            )

            if (uniform_color_scale) {

                x_lims <- global_x_lims
                y_lims <- global_y_lims

                x_breaks <- global_x_breaks
                y_breaks <- global_y_breaks

            } else {

                x_lims <- get_percentile_limits(
                    x = dat[[features[1]]],
                    bottom_percentile = x_bottom_percentile,
                    top_percentile = x_top_percentile,
                    feature_name = features[1]
                )

                y_lims <- get_percentile_limits(
                    x = dat[[features[2]]],
                    bottom_percentile = y_bottom_percentile,
                    top_percentile = y_top_percentile,
                    feature_name = features[2]
                )

                if (use_binned_colors) {
                    x_breaks <- make_bin_breaks(
                        limits = x_lims,
                        n_bins = n_bins
                    )

                    y_breaks <- make_bin_breaks(
                        limits = y_lims,
                        n_bins = n_bins
                    )
                } else {
                    x_breaks <- NULL
                    y_breaks <- NULL
                }
            }

            x_norm <- map_to_color_grid(
                x = dat[[features[1]]],
                limits = x_lims,
                side_length = side_length,
                n_bins = n_bins,
                bin_breaks = x_breaks
            )

            y_norm <- map_to_color_grid(
                x = dat[[features[2]]],
                limits = y_lims,
                side_length = side_length,
                n_bins = n_bins,
                bin_breaks = y_breaks
            )

            dat_norm <- cbind(x_norm, y_norm)

            feature_1_colors <- colorRampPalette(
                c(bottom_left, bottom_right)
            )(side_length)[x_norm]

            feature_2_colors <- colorRampPalette(
                c(bottom_left, top_left)
            )(side_length)[y_norm]

            blended_colors <- sapply(seq_len(nrow(dat_norm)), function(j) {
                col_grid[dat_norm[j, 1], dat_norm[j, 2]]
            })

            colors_list <- list(
                feature_1_colors,
                feature_2_colors,
                blended_colors
            )

            names(colors_list) <- c(
                features,
                paste0(features[1], "_", features[2])
            )

            # -----------------------------------------------------------------
            # Legend
            # -----------------------------------------------------------------

            if (use_binned_colors) {

                x_mids <- (head(x_breaks, -1) + tail(x_breaks, -1)) / 2
                y_mids <- (head(y_breaks, -1) + tail(y_breaks, -1)) / 2

                legend_grid <- expand.grid(
                    x_bin = seq_len(n_bins),
                    y_bin = seq_len(n_bins)
                )

                legend_grid[[features[1]]] <- x_mids[legend_grid$x_bin]
                legend_grid[[features[2]]] <- y_mids[legend_grid$y_bin]

                legend_grid$color <- sapply(seq_len(nrow(legend_grid)), function(j) {
                    col_grid[
                        legend_grid$x_bin[j],
                        legend_grid$y_bin[j]
                    ]
                })

                legend <- ggplot(
                    legend_grid,
                    aes(
                        x = .data[[features[1]]],
                        y = .data[[features[2]]],
                        fill = color
                    )
                ) +
                    geom_tile() +
                    scale_fill_identity() +
                    coord_cartesian(expand = FALSE) +
                    theme(
                        legend.position = "none",
                        aspect.ratio = 1,
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1)
                    ) +
                    xlab(ifelse(
                        is.null(feature_1_alt_name),
                        features[1],
                        feature_1_alt_name
                    )) +
                    ylab(ifelse(
                        is.null(feature_2_alt_name),
                        features[2],
                        feature_2_alt_name
                    ))

            } else {

                legend_grid <- expand.grid(
                    seq(
                        from = x_lims[1],
                        to = x_lims[2],
                        length.out = side_length
                    ),
                    seq(
                        from = y_lims[1],
                        to = y_lims[2],
                        length.out = side_length
                    )
                )

                colnames(legend_grid) <- features

                legend_grid$x_index <- rep(
                    seq_len(side_length),
                    times = side_length
                )

                legend_grid$y_index <- rep(
                    seq_len(side_length),
                    each = side_length
                )

                legend_grid$color <- sapply(seq_len(nrow(legend_grid)), function(j) {
                    col_grid[
                        legend_grid$x_index[j],
                        legend_grid$y_index[j]
                    ]
                })

                legend <- ggplot(
                    legend_grid,
                    aes(
                        x = .data[[features[1]]],
                        y = .data[[features[2]]],
                        color = color
                    )
                ) +
                    geom_point(shape = 15, size = 1.9) +
                    scale_color_identity() +
                    coord_cartesian(expand = FALSE) +
                    theme(
                        legend.position = "none",
                        aspect.ratio = 1,
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1)
                    ) +
                    xlab(ifelse(
                        is.null(feature_1_alt_name),
                        features[1],
                        feature_1_alt_name
                    )) +
                    ylab(ifelse(
                        is.null(feature_2_alt_name),
                        features[2],
                        feature_2_alt_name
                    ))
            }

        } else {

            if (top_right != "#FFFF00") {
                warning(paste(
                    "Cannot alter color in top right corner when",
                    "use_seurat_backend is TRUE."
                ))
            }

            vis_reduc <- cells_obj_sub@images[[i]]@coordinates[, c(3, 2)]
            colnames(vis_reduc) <- c("vis_1", "vis_2")
            vis_reduc$vis_2 <- -1 * vis_reduc$vis_2

            vis_reduc_mat <- as.matrix(vis_reduc)

            vis_reduc_obj <- CreateDimReducObject(
                embeddings = vis_reduc_mat,
                key = "vis_"
            )

            cells_obj_sub@reductions$vis <- vis_reduc_obj

            seurat_fp <- do.call(
                FeaturePlot,
                c(
                    list(
                        object = cells_obj_sub,
                        features = features,
                        reduction = "vis",
                        blend = TRUE,
                        cols = c(bottom_left, bottom_right, top_left),
                        combine = FALSE
                    ),
                    fp_extra_arguments
                )
            )

            colors_list <- lapply(seurat_fp[1:3], extract_colors_from_ggplot)

            names(colors_list) <- c(
                features,
                paste0(features[1], "_", features[2])
            )

            legend <- seurat_fp[[4]]
        }

        plot_list <- lapply(names(colors_list), function(x) {
            do.call(
                custom_color_SpatialDimPlot,
                c(
                    list(
                        cells_obj = cells_obj_sub,
                        image_name = i,
                        new_md_column_name = x,
                        colors_per_spot = colors_list[[x]]
                    ),
                    sfp_extra_arguments
                )
            )
        })

        plot_list[[4]] <- wrap_plots(
            ggplot() + theme_void(),
            legend,
            ggplot() + theme_void(),
            ncol = 1,
            heights = c(0.2, 0.6, 0.2)
        )

        plot_list_outer[[i]] <- plot_list
    }

    if (combine == FALSE) {
        return(plot_list_outer)
    } else {

        plot_list_outer <- lapply(plot_list_outer, function(p) {
            wrap_plots(
                p,
                nrow = 1,
                widths = c(0.28, 0.28, 0.28, 0.16)
            )
        })

        p <- wrap_plots(plot_list_outer, ncol = 1)

        return(p)
    }
}
