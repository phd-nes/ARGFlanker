#' Plot the gene neighbourhood of an ARG
#'
#' Draws the genes flanking one ARG along its contig, coloured by functional
#' category, as a quick visual check of a characterized neighbourhood.
#' Requires the \pkg{ggplot2} package.
#'
#' @param gr A classified [GenomicRanges::GRanges].
#' @param arg_id The `gene_id` of the ARG to centre on.
#' @param window Half-width in base pairs to show on each side (default
#'   `10000`).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- classifyGenes(importAnnotations(gff))
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'     plotNeighborhood(gr, arg_id = gr$gene_id[gr$category == "ARG"][1])
#' }
#' @export
plotNeighborhood <- function(gr, arg_id, window = 10000) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("plotNeighborhood() requires the 'ggplot2' package")

    gid <- as.character(S4Vectors::mcols(gr)$gene_id)
    i <- which(gid == arg_id)
    if (length(i) == 0) stop("arg_id '", arg_id, "' not found")
    i <- i[1]

    contig <- as.character(GenomeInfoDb::seqnames(gr))
    sel <- which(contig == contig[i] &
                 BiocGenerics::end(gr) >= BiocGenerics::start(gr[i]) - window &
                 BiocGenerics::start(gr) <= BiocGenerics::end(gr[i]) + window)

    df <- data.frame(
        gene_id = gid[sel],
        start = BiocGenerics::start(gr)[sel],
        end = BiocGenerics::end(gr)[sel],
        strand = as.character(BiocGenerics::strand(gr))[sel],
        category = as.character(S4Vectors::mcols(gr)$category)[sel],
        stringsAsFactors = FALSE)
    df$y <- 1
    df$is_focus <- df$gene_id == arg_id

    ggplot2::ggplot(df) +
        ggplot2::geom_rect(
            ggplot2::aes(xmin = start, xmax = end,
                         ymin = y - 0.3, ymax = y + 0.3,
                         fill = category),
            colour = "grey20") +
        ggplot2::geom_text(
            ggplot2::aes(x = (start + end) / 2, y = y + 0.45,
                         label = gene_id),
            size = 2.6, angle = 30, hjust = 0) +
        ggplot2::labs(
            title = paste0("Neighbourhood of ", arg_id,
                           " (", contig[i], ")"),
            x = "position (bp)", y = NULL, fill = "category") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_blank())
}
