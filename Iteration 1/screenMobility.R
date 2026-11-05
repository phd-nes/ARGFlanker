## MGE marker categories (everything that is an element marker, not the ARG
## itself and not background "other").
.mgeCategories <- function() {
    c("integrase", "transposase", "recombinase", "IS", "plasmid_rep")
}

#' Screen ARGs for co-located mobile genetic element markers
#'
#' For every ARG feature, finds MGE-marker genes lying on the *same contig*
#' within a distance window, and returns each ARG--MGE pair with the gap
#' between them. This is the contig-level ("physical linkage") evidence.
#'
#' Because short-read assembly often fragments an ARG away from its flanking
#' element, absence of a co-located marker here is weak evidence of
#' non-mobility; combine with [linkByReadPairs()] and [neighborhoodProfile()]
#' through [mobilityScore()].
#'
#' @param gr A classified [GenomicRanges::GRanges] (output of
#'   [classifyGenes()]) carrying a `category` column.
#' @param window Maximum gap, in base pairs, between an ARG and an MGE marker
#'   for them to count as linked (default `10000`).
#' @param mge_categories Character vector of categories treated as MGE markers.
#'
#' @return A [S4Vectors::DataFrame] with one row per ARG--MGE pair: `arg_id`,
#'   `arg_product`, `contig`, `mge_id`, `mge_category`, `mge_product`,
#'   `distance` and `source` (`"contig"`). Zero rows if nothing is linked.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' screenMobility(gr, window = 10000)
#' @export
screenMobility <- function(gr,
                           window = 10000,
                           mge_categories = .mgeCategories()) {
    if (!"category" %in% colnames(S4Vectors::mcols(gr)))
        stop("'gr' has no 'category' column; run classifyGenes() first")

    cat_vec <- as.character(S4Vectors::mcols(gr)$category)
    arg_idx <- which(cat_vec == "ARG")
    mge_idx <- which(cat_vec %in% mge_categories)

    empty <- S4Vectors::DataFrame(
        arg_id = character(0), arg_product = character(0),
        contig = character(0), mge_id = character(0),
        mge_category = character(0), mge_product = character(0),
        distance = integer(0), source = character(0))
    if (length(arg_idx) == 0 || length(mge_idx) == 0) return(empty)

    gid <- as.character(S4Vectors::mcols(gr)$gene_id)
    prod <- as.character(S4Vectors::mcols(gr)$product)
    contig <- as.character(GenomeInfoDb::seqnames(gr))

    rows <- vector("list", length(arg_idx))
    for (k in seq_along(arg_idx)) {
        i <- arg_idx[k]
        same <- mge_idx[contig[mge_idx] == contig[i]]
        if (length(same) == 0) next
        d <- IRanges::distance(gr[rep(i, length(same))], gr[same])
        keep <- which(!is.na(d) & d <= window)
        if (length(keep) == 0) next
        j <- same[keep]
        rows[[k]] <- S4Vectors::DataFrame(
            arg_id = gid[i], arg_product = prod[i], contig = contig[i],
            mge_id = gid[j], mge_category = cat_vec[j], mge_product = prod[j],
            distance = as.integer(d[keep]), source = "contig")
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(empty)
    do.call(rbind, rows)
}
