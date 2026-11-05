#' MobilityScreen class
#'
#' An S4 container holding the result of an ARGflanker mobility screen.
#'
#' @slot linkages A [S4Vectors::DataFrame-class] of ARG--MGE associations,
#'   one row per (ARG, MGE-marker) pair, with distance and evidence source.
#' @slot neighborhood A `data.frame` giving the enrichment of each gene
#'   category in the window around ARGs (the neighbourhood signature).
#' @slot scores A [S4Vectors::DataFrame-class] with one row per ARG and its
#'   combined mobility score.
#' @slot params A named `list` recording the parameters used for the screen.
#'
#' @name MobilityScreen-class
#' @rdname MobilityScreen-class
#' @exportClass MobilityScreen
setClass(
    "MobilityScreen",
    representation(
        linkages = "DataFrame",
        neighborhood = "data.frame",
        scores = "DataFrame",
        params = "list"
    ),
    prototype(
        linkages = S4Vectors::DataFrame(),
        neighborhood = data.frame(),
        scores = S4Vectors::DataFrame(),
        params = list()
    )
)

setValidity("MobilityScreen", function(object) {
    msg <- character(0)
    if (!is.list(object@params))
        msg <- c(msg, "'params' must be a list")
    if (length(msg)) msg else TRUE
})
