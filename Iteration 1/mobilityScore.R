## Bind linkage tables that may have different columns (contig vs readpair).
.bindLinkages <- function(...) {
    tabs <- list(...)
    tabs <- tabs[vapply(tabs, function(x) !is.null(x) && nrow(x) > 0, logical(1))]
    cols <- c("arg_id", "arg_product", "contig", "mge_id", "mge_category",
              "mge_product", "distance", "source", "support")
    if (length(tabs) == 0)
        return(S4Vectors::DataFrame(
            arg_id = character(0), arg_product = character(0),
            contig = character(0), mge_id = character(0),
            mge_category = character(0), mge_product = character(0),
            distance = integer(0), source = character(0),
            support = integer(0)))
    tabs <- lapply(tabs, function(tb) {
        if (!"support" %in% colnames(tb)) tb$support <- NA_integer_
        tb[, cols]
    })
    do.call(rbind, tabs)
}

#' Combine evidence into a per-ARG mobility score
#'
#' Integrates the three evidence sources into a single score in `[0, 1]` for
#' every ARG, using a soft-OR so that strong support from any one source raises
#' the score: contig proximity to an MGE marker, read-pair bridging across
#' contigs, and the neighbourhood signature from [neighborhoodProfile()].
#'
#' The score is a screening prior, not a probability of mobility. Confirm hits
#' with long reads, Hi-C, or targeted PCR before claiming an ARG is mobile.
#'
#' @param gr A classified [GenomicRanges::GRanges].
#' @param linkages A linkage [S4Vectors::DataFrame] from [screenMobility()]
#'   and/or [linkByReadPairs()] (combine with `rbind` is handled internally).
#' @param profile Neighbourhood signature from [neighborhoodProfile()].
#' @param window Neighbourhood half-width used for the signature (bp).
#' @param decay Distance decay constant for contig proximity (bp); smaller means
#'   proximity must be tighter to score high. Default `window / 3`.
#' @param support_scale Read-pair count at which bridging support is ~0.63 of
#'   its maximum (default `5`).
#' @param threshold Score at or above which an ARG is flagged `mobile`
#'   (default `0.5`).
#'
#' @return A [S4Vectors::DataFrame], one row per ARG, with component scores
#'   (`proximity_score`, `readpair_score`, `neighborhood_score`), the combined
#'   `score`, and a logical `mobile`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' lk <- screenMobility(gr)
#' prof <- neighborhoodProfile(gr)
#' mobilityScore(gr, lk, prof)
#' @export
mobilityScore <- function(gr, linkages, profile,
                          window = 10000,
                          decay = window / 3,
                          support_scale = 5,
                          threshold = 0.5) {
    cat_vec <- as.character(S4Vectors::mcols(gr)$category)
    gid <- as.character(S4Vectors::mcols(gr)$gene_id)
    prod <- as.character(S4Vectors::mcols(gr)$product)
    contig <- as.character(GenomeInfoDb::seqnames(gr))
    arg_idx <- which(cat_vec == "ARG")

    out <- S4Vectors::DataFrame(
        arg_id = gid[arg_idx], arg_product = prod[arg_idx],
        contig = contig[arg_idx])
    n <- length(arg_idx)
    out$n_contig_mge <- 0L
    out$nearest_distance <- NA_integer_
    out$n_readpair <- 0L
    out$readpair_support <- 0L
    out$proximity_score <- 0
    out$readpair_score <- 0
    out$neighborhood_score <- 0
    if (n == 0) {
        out$score <- numeric(0); out$mobile <- logical(0); return(out)
    }

    lk <- linkages
    if (is.null(lk) || nrow(lk) == 0) lk <- .bindLinkages()

    # Per-ARG neighbourhood score from the signature weights.
    wts <- setNames(profile$weight, profile$category)
    max_w <- sum(pmax(0, wts))
    args_gr <- gr[arg_idx]
    flanks <- GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(args_gr),
        ranges = IRanges::IRanges(
            start = pmax(1L, BiocGenerics::start(args_gr) - window),
            end = BiocGenerics::end(args_gr) + window),
        strand = "*")
    non_arg <- which(cat_vec != "ARG")

    for (k in seq_len(n)) {
        arg_id <- gid[arg_idx[k]]

        sub <- lk[!is.na(lk$arg_id) & lk$arg_id == arg_id, , drop = FALSE]
        contig_sub <- sub[sub$source == "contig", , drop = FALSE]
        rp_sub <- sub[sub$source == "readpair", , drop = FALSE]

        if (nrow(contig_sub) > 0) {
            out$n_contig_mge[k] <- length(unique(contig_sub$mge_id))
            nd <- suppressWarnings(min(contig_sub$distance, na.rm = TRUE))
            out$nearest_distance[k] <- as.integer(nd)
            out$proximity_score[k] <- exp(-nd / decay)
        }
        if (nrow(rp_sub) > 0) {
            out$n_readpair[k] <- length(unique(rp_sub$mge_id))
            sup <- sum(rp_sub$support, na.rm = TRUE)
            out$readpair_support[k] <- as.integer(sup)
            out$readpair_score[k] <- 1 - exp(-sup / support_scale)
        }

        if (length(non_arg) > 0 && max_w > 0) {
            ov <- IRanges::findOverlaps(gr[non_arg], flanks[k])
            near_cats <- unique(cat_vec[non_arg[S4Vectors::queryHits(ov)]])
            near_cats <- intersect(near_cats, names(wts))
            out$neighborhood_score[k] <- min(1, sum(wts[near_cats]) / max_w)
        }
    }

    # Soft-OR of the three components.
    p1 <- out$proximity_score
    p2 <- out$readpair_score
    p3 <- out$neighborhood_score
    out$score <- 1 - (1 - p1) * (1 - p2) * (1 - p3)
    out$mobile <- out$score >= threshold
    out
}

#' Run the full ARGflanker mobility screen
#'
#' Convenience constructor that runs [screenMobility()],
#' [neighborhoodProfile()], optional [linkByReadPairs()] and [mobilityScore()],
#' and returns a [MobilityScreen-class] object.
#'
#' @param gr A classified [GenomicRanges::GRanges] (see [classifyGenes()]).
#' @param pairs Optional read-pair `data.frame` for [linkByReadPairs()].
#' @param window Distance window / neighbourhood half-width in bp.
#' @param threshold Mobility-score cutoff for the `mobile` flag.
#' @param ... Further arguments passed to [mobilityScore()].
#'
#' @return A [MobilityScreen-class] object.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' res <- MobilityScreen(gr)
#' res
#' @export
MobilityScreen <- function(gr, pairs = NULL, window = 10000,
                           threshold = 0.5, ...) {
    contig_lk <- screenMobility(gr, window = window)
    rp_lk <- if (!is.null(pairs)) linkByReadPairs(pairs) else NULL
    lk <- .bindLinkages(contig_lk, rp_lk)
    prof <- neighborhoodProfile(gr, window = window)
    sc <- mobilityScore(gr, lk, prof, window = window,
                        threshold = threshold, ...)
    methods::new("MobilityScreen",
        linkages = lk, neighborhood = prof, scores = sc,
        params = list(window = window, threshold = threshold,
                      used_readpairs = !is.null(pairs)))
}
