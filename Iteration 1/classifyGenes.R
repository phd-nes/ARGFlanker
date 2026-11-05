#' Default keyword patterns for gene categories
#'
#' Returns the regular-expression patterns used by [classifyGenes()] to assign
#' each predicted gene to a functional category from its product/annotation
#' text. The patterns are deliberately permissive and are meant to be edited:
#' in a real pipeline you would replace keyword matching with hits against a
#' curated reference (e.g. CARD/Resfams for ARGs, ISfinder / mobileOG for MGE
#' markers).
#'
#' @return A named `list` of character vectors. Names are category labels; the
#'   order defines precedence (earlier categories win when a product matches
#'   more than one pattern). The final category, `"other"`, matches nothing and
#'   acts as the fallback.
#'
#' @examples
#' defaultPatterns()
#' @export
defaultPatterns <- function() {
    list(
        ARG = c("resistance", "beta-?lactamase", "\\bbla[A-Z]", "aminoglycoside",
                "tetracycline", "\\btet[A-Z]", "macrolide", "\\berm[A-Z]",
                "chloramphenicol", "\\bcat\\b", "sulfonamide", "\\bsul[0-9]",
                "vancomycin", "\\bvan[A-Z]", "efflux.*antibiotic", "\\bqnr",
                "\\bmcr-?[0-9]", "\\bdfr[A-Z]", "\\baac\\b", "\\baph\\b"),
        integrase = c("integrase", "\\bint[A-Z]?[0-9]?\\b", "integron"),
        transposase = c("transposase", "\\btnp", "\\btn[0-9]", "Tn3"),
        recombinase = c("recombinase", "resolvase", "\\bxer[CD]"),
        IS = c("insertion sequence", "\\bIS[0-9]", "\\bISCR"),
        plasmid_rep = c("plasmid", "\\brep[AB]?\\b", "replication initiat",
                        "mobilization", "\\bmob[A-Z]", "relaxase", "\\btra[A-Z]",
                        "conjugal", "conjugative", "type IV secretion"),
        other = character(0)
    )
}

#' Classify predicted genes into functional categories
#'
#' Assigns each feature to one of the ARGflanker categories (ARG, integrase,
#' transposase, recombinase, IS, plasmid_rep, other) from a free-text
#' annotation column, using keyword patterns.
#'
#' @param gr A [GenomicRanges::GRanges] (e.g. from [importAnnotations()]) or a
#'   `data.frame` with an annotation column.
#' @param annotation_col Name of the column holding the product/annotation text.
#' @param patterns A named list of regex patterns; see [defaultPatterns()].
#'   Categories are applied in list order and the first match wins.
#' @param ignore_case Logical; case-insensitive matching (default `TRUE`).
#'
#' @return The input with an added/overwritten `category` factor column (for a
#'   `GRanges`, added to `mcols`).
#'
#' @examples
#' gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
#' gr <- importAnnotations(gff)
#' gr <- classifyGenes(gr)
#' table(gr$category)
#' @export
classifyGenes <- function(gr,
                          annotation_col = "product",
                          patterns = defaultPatterns(),
                          ignore_case = TRUE) {
    cats <- names(patterns)
    if (is(gr, "GRanges")) {
        md <- S4Vectors::mcols(gr)
        if (!annotation_col %in% colnames(md))
            stop("annotation column '", annotation_col,
                 "' not found in mcols(gr)")
        text <- as.character(md[[annotation_col]])
    } else if (is.data.frame(gr)) {
        if (!annotation_col %in% colnames(gr))
            stop("annotation column '", annotation_col, "' not found")
        text <- as.character(gr[[annotation_col]])
    } else {
        stop("'gr' must be a GRanges or data.frame")
    }
    text[is.na(text)] <- ""

    category <- rep("other", length(text))
    assigned <- logical(length(text))
    for (cat in cats) {
        pats <- patterns[[cat]]
        if (length(pats) == 0) next
        rx <- paste(pats, collapse = "|")
        hit <- grepl(rx, text, ignore.case = ignore_case, perl = TRUE)
        take <- hit & !assigned
        category[take] <- cat
        assigned <- assigned | take
    }
    category <- factor(category, levels = cats)

    if (is(gr, "GRanges")) {
        S4Vectors::mcols(gr)$category <- category
    } else {
        gr$category <- category
    }
    gr
}
