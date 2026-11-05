#' ARGflanker: characterizing the gene neighbourhoods of antibiotic
#' resistance genes
#'
#' ARGflanker describes and scores the genes that flank antibiotic resistance
#' genes (ARGs) in short-read metagenomes. Its primary output is a
#' characterization of each ARG's neighbourhood - the categories, distances and
#' orientations of nearby genes - from which it derives a neighbourhood
#' signature and a per-ARG mobility score that also incorporates cross-contig
#' read-pair evidence.
#'
#' Core functions:
#' \describe{
#'   \item{[importAnnotations()]}{read ORF predictions (GFF3 or data.frame).}
#'   \item{[classifyGenes()]}{assign genes to functional categories.}
#'   \item{[characterizeNeighbors()]}{per-ARG table of flanking genes.}
#'   \item{[summarizeNeighbors()]}{one-row-per-ARG neighbourhood summary.}
#'   \item{[neighborhoodProfile()]}{dataset-wide enrichment signature.}
#'   \item{[screenMobility()] / [linkByReadPairs()]}{ARG-MGE linkage evidence.}
#'   \item{[mobilityScore()] / [MobilityScreen()]}{combined scoring.}
#' }
#'
#' @keywords internal
"_PACKAGE"
