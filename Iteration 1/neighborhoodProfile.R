#' Learn the gene-category signature around ARGs
#'
#' Builds the "guilt by association" signature: for a window around every ARG,
#' counts the functional categories of neighbouring genes and tests whether each
#' category is enriched in ARG neighbourhoods relative to its genome-wide
#' frequency (one-sided hypergeometric test). The resulting weights feed
#' [mobilityScore()] so that an ARG surrounded by the categories that *usually*
#' flank resistance genes scores higher, even without a single confirmed MGE
#' hit.
#'
#' @param gr A classified [GenomicRanges::GRanges] (see [classifyGenes()]).
#' @param window Half-width, in base pairs, of the neighbourhood searched on
#'   each side of an ARG (default `10000`).
#' @param pseudocount Added to counts to stabilise ratios (default `1`).
#'
#' @return A `data.frame` with one row per non-ARG category: `category`,
#'   `n_near_arg` (occurrences within windows of ARGs), `n_total`,
#'   `frac_near`, `frac_background`, `enrichment` (frac_near / frac_background)
#'   and `p_value`. A `weight` column (clipped log2 enrichment, >= 0) is added
#'   for scoring.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' neighborhoodProfile(gr)
#' @export
neighborhoodProfile <- function(gr, window = 10000, pseudocount = 1) {
    if (!"category" %in% colnames(S4Vectors::mcols(gr)))
        stop("'gr' has no 'category' column; run classifyGenes() first")

    cat_vec <- as.character(S4Vectors::mcols(gr)$category)
    all_cats <- setdiff(unique(cat_vec), "ARG")
    arg_idx <- which(cat_vec == "ARG")

    # Expand each ARG into a flanking window and find neighbouring genes.
    near_counts <- setNames(rep(0L, length(all_cats)), all_cats)
    if (length(arg_idx) > 0) {
        args_gr <- gr[arg_idx]
        flanks <- GenomicRanges::GRanges(
            seqnames = GenomeInfoDb::seqnames(args_gr),
            ranges = IRanges::IRanges(
                start = pmax(1L, BiocGenerics::start(args_gr) - window),
                end = BiocGenerics::end(args_gr) + window),
            strand = "*")
        non_arg <- which(cat_vec != "ARG")
        if (length(non_arg) > 0) {
            hits <- IRanges::findOverlaps(gr[non_arg], flanks)
            near <- non_arg[S4Vectors::queryHits(hits)]
            near <- unique(near)  # count each neighbour once
            tt <- table(factor(cat_vec[near], levels = all_cats))
            near_counts[names(tt)] <- as.integer(tt)
        }
    }

    total_counts <- table(factor(cat_vec, levels = c("ARG", all_cats)))
    n_total_genes <- length(cat_vec)
    n_near_total <- sum(near_counts)

    frac_near <- (near_counts + pseudocount) /
        (n_near_total + pseudocount * length(all_cats))
    frac_bg <- (as.integer(total_counts[all_cats]) + pseudocount) /
        (n_total_genes + pseudocount * length(all_cats))
    enrichment <- frac_near / frac_bg

    # Hypergeometric: P(X >= observed) of drawing this category into the
    # neighbourhood sample from the whole gene pool.
    pval <- vapply(all_cats, function(ct) {
        white <- as.integer(total_counts[ct])
        black <- n_total_genes - white
        drawn <- n_near_total
        obs <- near_counts[[ct]]
        if (drawn == 0 || white == 0) return(1)
        stats::phyper(obs - 1L, white, black, drawn, lower.tail = FALSE)
    }, numeric(1))

    weight <- pmax(0, log2(enrichment))

    data.frame(
        category = all_cats,
        n_near_arg = as.integer(near_counts[all_cats]),
        n_total = as.integer(total_counts[all_cats]),
        frac_near = as.numeric(frac_near),
        frac_background = as.numeric(frac_bg),
        enrichment = as.numeric(enrichment),
        p_value = as.numeric(pval),
        weight = as.numeric(weight),
        row.names = NULL,
        stringsAsFactors = FALSE)
}
