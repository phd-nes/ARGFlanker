#' Import gene predictions and annotations
#'
#' Reads predicted genes (ORFs) with their contig coordinates into a
#' [GenomicRanges::GRanges]. Accepts a GFF3 file (as produced by Prodigal +
#' a functional annotator, or by prokka/bakta) or a `data.frame`.
#'
#' Short-read workflows typically assemble reads into contigs, call ORFs on the
#' contigs, and annotate them; ARGflanker consumes that per-contig ORF table.
#' Each `seqnames` value is a contig and the ranges are ORF positions along it.
#'
#' @param x Path to a GFF3 file, or a `data.frame` with columns `contig`,
#'   `start`, `end`, `strand`, `gene_id` and an annotation column.
#' @param annotation_col For a `data.frame`, the name of the product/annotation
#'   column (default `"product"`). For GFF3, the attribute searched for product
#'   text (default tries `product`, then `Name`, then `ID`).
#' @param feature_type For GFF3, keep only rows with this `type` (default
#'   `"CDS"`); set `NULL` to keep all.
#'
#' @return A `GRanges` with one range per gene and a metadata column holding the
#'   annotation text (named `product`).
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- importAnnotations(gff)
#' gr
#' @export
importAnnotations <- function(x,
                              annotation_col = "product",
                              feature_type = "CDS") {
    if (is.character(x) && length(x) == 1L) {
        gr <- rtracklayer::import(x)
        if (!is.null(feature_type) && "type" %in% colnames(S4Vectors::mcols(gr)))
            gr <- gr[as.character(gr$type) %in% feature_type]
        md <- S4Vectors::mcols(gr)
        prod <- NULL
        for (cand in c(annotation_col, "product", "Name", "Note", "ID")) {
            if (cand %in% colnames(md)) {
                prod <- as.character(md[[cand]])
                break
            }
        }
        if (is.null(prod)) prod <- rep(NA_character_, length(gr))
        gid <- if ("ID" %in% colnames(md)) as.character(md$ID)
               else if ("gene_id" %in% colnames(md)) as.character(md$gene_id)
               else paste0("gene", seq_along(gr))
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            gene_id = gid, product = prod)
        return(gr)
    }

    if (is.data.frame(x)) {
        req <- c("contig", "start", "end")
        miss <- setdiff(req, colnames(x))
        if (length(miss))
            stop("data.frame is missing columns: ", paste(miss, collapse = ", "))
        str <- if ("strand" %in% colnames(x)) as.character(x$strand) else "*"
        str[is.na(str) | !str %in% c("+", "-", "*")] <- "*"
        gid <- if ("gene_id" %in% colnames(x)) as.character(x$gene_id)
               else paste0("gene", seq_len(nrow(x)))
        prod <- if (annotation_col %in% colnames(x))
            as.character(x[[annotation_col]]) else NA_character_
        gr <- GenomicRanges::GRanges(
            seqnames = as.character(x$contig),
            ranges = IRanges::IRanges(start = as.integer(x$start),
                                      end = as.integer(x$end)),
            strand = str)
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            gene_id = gid, product = prod)
        return(gr)
    }

    stop("'x' must be a GFF3 file path or a data.frame")
}
