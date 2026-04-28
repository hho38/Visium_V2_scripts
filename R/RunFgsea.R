#' Run preranked FGSEA from a differential expression table
#'
#' Performs fast gene set enrichment analysis (FGSEA) using a differential
#' expression results table as input. Genes are ranked using a signed
#' transformation of p-values:
#'
#' \deqn{sign(logFC) * -log10(pvalue)}
#'
#' This ranking encodes direction from the log fold change and significance
#' from the p-value. If pathway gene sets are not supplied, the function
#' retrieves MSigDB gene sets using \pkg{msigdbr}.
#'
#' Duplicate genes are resolved by keeping the entry with the largest absolute
#' ranking statistic.
#'
#' @param de A data frame containing differential expression results.
#' @param pathways An optional named list of pathways, where each element is a
#'   character vector of gene symbols. If `NULL`, pathways are retrieved from
#'   \pkg{msigdbr} using `species`, `collection`, and `subcollection`.
#' @param species Species name passed to [msigdbr::msigdbr()], such as
#'   `"Mus musculus"` or `"Homo sapiens"`.
#' @param collection MSigDB collection passed to [msigdbr::msigdbr()].
#'   Common examples include `"H"` for Hallmark and `"C5"` for Gene Ontology.
#' @param subcollection Optional MSigDB subcollection passed to
#'   [msigdbr::msigdbr()]. Default is `NULL`.
#' @param gene_col Optional column name containing gene symbols. If `NULL`,
#'   row names of `de` are used.
#' @param lfc_col Optional column name containing log fold change values. If
#'   `NULL`, the function looks for `"avg_log2FC"` first, then `"avg_logFC"`.
#' @param pval_col Column name containing p-values. Default is `"p_val"`.
#' @param minSize Minimum pathway size passed to [fgsea::fgsea()].
#' @param maxSize Maximum pathway size passed to [fgsea::fgsea()].
#' @param jitter Numeric value controlling the spacing added to ranks to break
#'   ties. Default is `0`, which leaves ranks unchanged apart from the initial
#'   `1e-12` offset used in the function.
#' @param nPermSimple Number of simple permutations used by [fgsea::fgsea()].
#'   Increasing this can improve p-value estimation at the cost of runtime. Default is NA to call fgseaMultilevel.
#' @param plot Logical; if `TRUE`, plots the ranked statistic distribution
#'   before running FGSEA.
#'
#' @return A data frame of FGSEA results sorted first by adjusted p-value
#'   (`padj`) and then by absolute normalized enrichment score (`NES`).
#'
#' @details
#' Rows with missing gene names, missing fold changes, missing p-values, or
#' non-finite rank values are removed before enrichment analysis. P-values less
#' than or equal to zero are replaced with `.Machine$double.xmin` to avoid
#' infinite values when computing `-log10(pvalue)`.
#'
#' The function sets `set.seed(123)` internally to improve reproducibility of
#' FGSEA results.
#'
#' Note that the default ranking statistic,
#' `sign(logFC) * -log10(pvalue)`, does not use the magnitude of the log fold
#' change. Alternative ranking strategies, such as using `logFC` alone or
#' `logFC * -log10(pvalue)`, may be more appropriate depending on the analysis.
#'
#' @examples
#' \dontrun{
#' fgsea_res <- run_fgsea_from_de(
#'   de = osteo_lung2_v_lung1_df,
#'   species = "Mus musculus",
#'   collection = "H",
#'   subcollection = NULL,
#'   jitter = 1e-12
#' )
#'
#' head(fgsea_res)
#' }
#'
#' @export
run_fgsea_from_de <- function(
  de,
  pathways = NULL,
  species = "Mus musculus",
  collection = "H",
  subcollection = NULL,
  gene_col = NULL,
  lfc_col = NULL,
  pval_col = "p_val",
  minSize = 10,
  maxSize = 500, 
  jitter = 0,
  nPermSimple = NA,
  plot = FALSE
) {
  if (!is.data.frame(de)) stop("'de' must be a data.frame")
  set.seed(123)

  if (is.null(lfc_col)) {
    if ("avg_log2FC" %in% colnames(de)) {
      lfc_col <- "avg_log2FC"
    } else if ("avg_logFC" %in% colnames(de)) {
      lfc_col <- "avg_logFC"
    } else {
      stop("Could not find logFC column. Expected 'avg_log2FC' or 'avg_logFC'.")
    }
  }

  if (!lfc_col %in% colnames(de)) stop(sprintf("Column '%s' not found", lfc_col))
  if (!pval_col %in% colnames(de)) stop(sprintf("Column '%s' not found", pval_col))

  genes <- if (is.null(gene_col)) rownames(de) else de[[gene_col]]
  if (is.null(genes)) stop("Gene names not found. Use NULL for rownames(de) or provide gene_col.")

  if (is.null(pathways)) {
    msig_df <- msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
    pathways <- split(msig_df$gene_symbol, msig_df$gs_name)
  }

  df <- data.frame(
    gene = as.character(genes),
    lfc = de[[lfc_col]],
    pval = de[[pval_col]],
    stringsAsFactors = FALSE
  )

  # Keep rows with valid gene names.
  df <- df[!is.na(df$gene) & df$gene != "", , drop = FALSE]

  # Keep rows with non-missing fold change and p-value.
  df <- df[!is.na(df$lfc) & !is.na(df$pval), , drop = FALSE]

  # Replace zero or negative p-values with the smallest positive double to
  # avoid infinite values when computing -log10(p-value).
  df$pval[df$pval <= 0] <- .Machine$double.xmin

  # Compute the preranked statistic using the sign of the fold change and the
  # significance of the p-value.
  #df$rank <- sign(df$lfc) * -log10(df$pval)
  #df$rank <- df$lfc
  #df$rank <- df$lfc * -log10(pmax(df$pval, .Machine$double.xmin))
  df$rank <- sign(df$lfc) * -log10(pmax(df$pval, .Machine$double.xmin))

  # Remove rows with non-finite rank values.
  df <- df[is.finite(df$rank), , drop = FALSE]
  if (nrow(df) == 0) stop("No valid genes remain after filtering")

  # Resolve duplicate genes by keeping the row with the largest absolute rank.
  df <- df[order(abs(df$rank), decreasing = TRUE), ]
  df <- df[!duplicated(df$gene), ]

  # Create a named numeric vector required by fgsea, with gene symbols as names
  # and ranking statistics as values.
  ranks <- df$rank
  
  # adds a tiny unique increment,or jitter, to every gene in the list. Ex 1e-12
  #ranks <- ranks + seq(1e-12, length.out = length(ranks), by = jitter)
  # Add a tiny deterministic offset only within tied rank groups.
  # This breaks ties without changing ranks that are already unique.
  if (jitter > 0) {
    ranks <- ave(
      ranks,
      ranks,
      FUN = function(x) {
        if (length(x) == 1) {
          x
        } else {
          x + seq(0, by = jitter, length.out = length(x))
        }
      }
    )
  }
  
  names(ranks) <- df$gene

  # Sort ranks from highest to lowest for preranked GSEA.
  ranks <- sort(ranks, decreasing = TRUE)
  
  if (plot == TRUE) {
    plot(ranks)
  }
  
  if (!is.na(nPermSimple)){
    #This will run fgseaSimple. Useful when pathways of interest have NA values
    fgsea_res <- fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = minSize,
      maxSize = maxSize, 
      nPermSimple = nPermSimple
    )
  }else{
    #This will run the default modern fgseaMultilevel
      fgsea_res <- fgsea::fgsea(
      pathways = pathways,
      stats = ranks,
      minSize = minSize,
      maxSize = maxSize,
    )
  }

  # Return results ordered by adjusted p-value, then by enrichment magnitude.
  fgsea_res[order(fgsea_res$padj, -abs(fgsea_res$NES)), ]
}