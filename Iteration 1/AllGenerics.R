#' Accessors for MobilityScreen objects
#'
#' @param object A [MobilityScreen-class] object.
#' @return `linkages()` returns the ARG--MGE linkage table; `neighborhood()`
#'   returns the neighbourhood-signature data frame; `screenParams()` returns
#'   the parameter list used for the screen.
#' @name MobilityScreen-accessors
#' @rdname MobilityScreen-accessors
#' @export
setGeneric("linkages", function(object) standardGeneric("linkages"))

#' @rdname MobilityScreen-accessors
#' @export
setGeneric("neighborhood", function(object) standardGeneric("neighborhood"))

#' @rdname MobilityScreen-accessors
#' @export
setGeneric("screenParams", function(object) standardGeneric("screenParams"))

#' @rdname MobilityScreen-accessors
setMethod("linkages", "MobilityScreen", function(object) object@linkages)

#' @rdname MobilityScreen-accessors
setMethod("neighborhood", "MobilityScreen", function(object) object@neighborhood)

#' @rdname MobilityScreen-accessors
setMethod("screenParams", "MobilityScreen", function(object) object@params)

#' @rdname MobilityScreen-class
#' @param object A [MobilityScreen-class] object.
#' @importFrom methods show
setMethod("show", "MobilityScreen", function(object) {
    cat("MobilityScreen\n")
    sc <- object@scores
    narg <- nrow(sc)
    nmob <- if (narg > 0 && "mobile" %in% colnames(sc))
        sum(sc$mobile, na.rm = TRUE) else NA_integer_
    cat("  ARGs screened:      ", narg, "\n", sep = "")
    cat("  flagged mobile:     ", nmob, "\n", sep = "")
    cat("  ARG-MGE linkages:   ", nrow(object@linkages), "\n", sep = "")
    win <- object@params[["window"]]
    if (!is.null(win)) cat("  window (bp):        ", win, "\n", sep = "")
    invisible(object)
})
