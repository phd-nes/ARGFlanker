#' Link ARGs to MGE markers across contigs using read pairs
#'
#' When an ARG and an MGE marker fall on *separate* contigs, assembly has lost
#' their physical connection, but paired-end reads whose two mates map to the
#' two different genes still testify that the genes were close on the original
#' DNA molecule. This function tallies such bridging pairs into ARG--MGE
#' linkages, recovering evidence that contig-level screening misses.
#'
#' @param pairs A `data.frame` of mapped read pairs with, for each mate, the
#'   gene it hit and that gene's category. Required columns: `gene_a`,
#'   `category_a`, `gene_b`, `category_b`.
#' @param min_pairs Minimum number of bridging read pairs required to report a
#'   linkage (default `2`).
#' @param mge_categories Categories treated as MGE markers.
#'
#' @return A [S4Vectors::DataFrame] of ARG--MGE linkages with `source =
#'   "readpair"` and a `support` column giving the bridging-pair count.
#'   Coordinate columns (`contig`, `distance`) are `NA` because the link spans
#'   contigs.
#'
#' @examples
#' pairs <- data.frame(
#'     gene_a = c("argA", "argA", "argA", "argB"),
#'     category_a = c("ARG", "ARG", "ARG", "ARG"),
#'     gene_b = c("tnpX", "tnpX", "intY", "housekeep1"),
#'     category_b = c("transposase", "transposase", "integrase", "other"))
#' linkByReadPairs(pairs, min_pairs = 1)
#' @export
linkByReadPairs <- function(pairs,
                            min_pairs = 2,
                            mge_categories = .mgeCategories()) {
    req <- c("gene_a", "category_a", "gene_b", "category_b")
    miss <- setdiff(req, colnames(pairs))
    if (length(miss))
        stop("'pairs' is missing columns: ", paste(miss, collapse = ", "))

    ca <- as.character(pairs$category_a)
    cb <- as.character(pairs$category_b)

    # Orient each pair so the ARG is side 'a' and the MGE marker is side 'b'.
    a_is_arg <- ca == "ARG" & cb %in% mge_categories
    b_is_arg <- cb == "ARG" & ca %in% mge_categories

    arg_gene <- c(as.character(pairs$gene_a)[a_is_arg],
                  as.character(pairs$gene_b)[b_is_arg])
    mge_gene <- c(as.character(pairs$gene_b)[a_is_arg],
                  as.character(pairs$gene_a)[b_is_arg])
    mge_cat  <- c(cb[a_is_arg], ca[b_is_arg])

    empty <- S4Vectors::DataFrame(
        arg_id = character(0), arg_product = character(0),
        contig = character(0), mge_id = character(0),
        mge_category = character(0), mge_product = character(0),
        distance = integer(0), source = character(0), support = integer(0))
    if (length(arg_gene) == 0) return(empty)

    key <- paste(arg_gene, mge_gene, mge_cat, sep = "\r")
    tab <- table(key)
    tab <- tab[tab >= min_pairs]
    if (length(tab) == 0) return(empty)

    parts <- do.call(rbind, strsplit(names(tab), "\r", fixed = TRUE))
    S4Vectors::DataFrame(
        arg_id = parts[, 1],
        arg_product = NA_character_,
        contig = NA_character_,
        mge_id = parts[, 2],
        mge_category = parts[, 3],
        mge_product = NA_character_,
        distance = NA_integer_,
        source = "readpair",
        support = as.integer(tab))
}
