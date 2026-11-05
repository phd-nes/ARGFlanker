#' Characterize the genes neighbouring each ARG
#'
#' The core ARGflanker operation: for every ARG, describe the genes that flank
#' it on the same contig within a window. For each neighbour it records the
#' functional category, the gap to the ARG, whether it lies upstream or
#' downstream, the rank (1 = nearest), and whether it shares the ARG's strand.
#' This is the raw "what lives next to resistance genes" table that every other
#' ARGflanker function is built on.
#'
#' Upstream/downstream are reported relative to the ARG's own strand when the
#' ARG is on `+`/`-`; for unstranded ARGs they fall back to genomic coordinate
#' order (lower start = "upstream").
#'
#' @param gr A classified [GenomicRanges::GRanges] (see [classifyGenes()]) with
#'   `gene_id`, `product` and `category` columns.
#' @param window Half-width in base pairs searched on each side of the ARG
#'   (default `10000`).
#' @param include_other Keep background (`other`) neighbours (default `TRUE`).
#'   Set `FALSE` to restrict to ARG/MGE categories.
#' @param max_rank Optional cap on the number of neighbours reported per side
#'   (nearest first); `NULL` keeps all within the window.
#'
#' @return A [S4Vectors::DataFrame], one row per (ARG, neighbour) pair:
#'   `arg_id`, `arg_product`, `arg_class` (if present), `contig`,
#'   `neighbor_id`, `neighbor_product`, `neighbor_category`, `distance`,
#'   `side` (`"upstream"`/`"downstream"`), `rank`, `same_strand`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' nb <- characterizeNeighbors(gr, window = 10000)
#' head(nb)
#' @export
characterizeNeighbors <- function(gr, window = 10000,
                                  include_other = TRUE,
                                  max_rank = NULL) {
    if (!"category" %in% colnames(S4Vectors::mcols(gr)))
        stop("'gr' has no 'category' column; run classifyGenes() first")

    cat_vec <- as.character(S4Vectors::mcols(gr)$category)
    gid <- as.character(S4Vectors::mcols(gr)$gene_id)
    prod <- as.character(S4Vectors::mcols(gr)$product)
    contig <- as.character(GenomeInfoDb::seqnames(gr))
    starts <- BiocGenerics::start(gr)
    ends <- BiocGenerics::end(gr)
    strands <- as.character(BiocGenerics::strand(gr))
    arg_class <- if ("arg_class" %in% colnames(S4Vectors::mcols(gr)))
        as.character(S4Vectors::mcols(gr)$arg_class) else NA_character_

    arg_idx <- which(cat_vec == "ARG")
    empty <- S4Vectors::DataFrame(
        arg_id = character(0), arg_product = character(0),
        arg_class = character(0), contig = character(0),
        neighbor_id = character(0), neighbor_product = character(0),
        neighbor_category = character(0), distance = integer(0),
        side = character(0), rank = integer(0), same_strand = logical(0))
    if (length(arg_idx) == 0) return(empty)

    rows <- vector("list", length(arg_idx))
    for (k in seq_along(arg_idx)) {
        i <- arg_idx[k]
        cand <- which(contig == contig[i] & seq_along(gr) != i)
        if (!include_other) cand <- cand[cat_vec[cand] != "other"]
        if (length(cand) == 0) next

        d <- IRanges::distance(gr[rep(i, length(cand))], gr[cand])
        keep <- which(!is.na(d) & d <= window)
        if (length(keep) == 0) next
        cand <- cand[keep]; d <- d[keep]

        # Side relative to the ARG, oriented by ARG strand where defined.
        downstream_coord <- starts[cand] >= ends[i]
        if (strands[i] == "-") {
            side <- ifelse(downstream_coord, "upstream", "downstream")
        } else {
            side <- ifelse(downstream_coord, "downstream", "upstream")
        }
        same_strand <- strands[cand] == strands[i] & strands[i] %in% c("+", "-")

        df <- S4Vectors::DataFrame(
            arg_id = gid[i], arg_product = prod[i],
            arg_class = arg_class[i], contig = contig[i],
            neighbor_id = gid[cand], neighbor_product = prod[cand],
            neighbor_category = cat_vec[cand],
            distance = as.integer(d), side = side,
            rank = NA_integer_, same_strand = same_strand)

        # Rank per side, nearest = 1.
        ord <- order(df$side, df$distance)
        df <- df[ord, , drop = FALSE]
        df$rank <- ave_rank(df$side, df$distance)
        if (!is.null(max_rank)) df <- df[df$rank <= max_rank, , drop = FALSE]
        rows[[k]] <- df
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(empty)
    do.call(rbind, rows)
}

## within-side rank by distance
ave_rank <- function(side, distance) {
    r <- integer(length(side))
    for (s in unique(side)) {
        idx <- which(side == s)
        r[idx] <- rank(distance[idx], ties.method = "first")
    }
    as.integer(r)
}

#' Summarize ARG neighbourhoods
#'
#' Collapses the per-neighbour table from [characterizeNeighbors()] into one row
#' per ARG: counts of each neighbour category, the nearest MGE marker and its
#' distance, and a compact description of the local gene context.
#'
#' @param neighbors A [S4Vectors::DataFrame] from [characterizeNeighbors()].
#' @param mge_categories Categories treated as MGE markers.
#'
#' @return A [S4Vectors::DataFrame], one row per ARG, with `n_neighbors`, a
#'   count column per category, `nearest_mge_category`, `nearest_mge_distance`
#'   and `has_mge`.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' summarizeNeighbors(characterizeNeighbors(gr))
#' @export
summarizeNeighbors <- function(neighbors, mge_categories = .mgeCategories()) {
    if (nrow(neighbors) == 0)
        return(S4Vectors::DataFrame(
            arg_id = character(0), n_neighbors = integer(0),
            nearest_mge_category = character(0),
            nearest_mge_distance = integer(0), has_mge = logical(0)))

    args <- unique(neighbors$arg_id)
    all_cats <- sort(unique(as.character(neighbors$neighbor_category)))

    res <- lapply(args, function(a) {
        sub <- neighbors[neighbors$arg_id == a, , drop = FALSE]
        counts <- as.list(table(factor(as.character(sub$neighbor_category),
                                       levels = all_cats)))
        mge <- sub[sub$neighbor_category %in% mge_categories, , drop = FALSE]
        if (nrow(mge) > 0) {
            j <- which.min(mge$distance)
            nmc <- as.character(mge$neighbor_category[j])
            nmd <- as.integer(mge$distance[j])
        } else {
            nmc <- NA_character_; nmd <- NA_integer_
        }
        base <- S4Vectors::DataFrame(
            arg_id = a, n_neighbors = nrow(sub),
            nearest_mge_category = nmc, nearest_mge_distance = nmd,
            has_mge = nrow(mge) > 0)
        for (ct in all_cats)
            base[[paste0("n_", ct)]] <- as.integer(counts[[ct]])
        base
    })
    do.call(rbind, res)
}
