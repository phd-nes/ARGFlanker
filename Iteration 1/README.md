# ARGflanker

**Characterizing the gene neighbourhoods of antibiotic resistance genes in
short-read metagenomes.**

ARGflanker describes what sits next to each antibiotic resistance gene (ARG) in
a metagenome — the functional categories, distances and orientations of the
flanking genes — and uses that context to flag ARGs that are likely to be
mobile. It is aimed at environmental resistomes (e.g. soil and airborne/dust
samples) where ARGs are often divergent and where short-read assembly breaks the
physical link between an ARG and its mobile genetic element (MGE).

> **Scope.** This is a hypothesis-generating screen. A high mobility score marks
> an ARG for experimental follow-up (long-read sequencing, Hi-C, targeted PCR,
> functional assays) — it is not proof of mobility or of conferred resistance.

## Why neighbourhoods

An ARG on its own tells you resistance *could* be present. Whether it can
*spread* depends on its context: integrases, transposases, recombinases,
insertion sequences and plasmid/conjugation genes nearby are the signature of a
mobilizable ARG. Because co-assembly frequently separates an ARG from these
markers, ARGflanker does not rely on a single contig joining them. It combines:

1. **Per-ARG neighbourhood characterization** — `characterizeNeighbors()`
2. **A dataset-wide enrichment signature** — `neighborhoodProfile()`
3. **On-contig linkage** — `screenMobility()`
4. **Cross-contig read-pair bridging** — `linkByReadPairs()`
5. **A combined, per-ARG mobility score** — `mobilityScore()` / `MobilityScreen()`

## Installation

```r
# install.packages("BiocManager")
# BiocManager::install("ARGflanker")      # once accepted to Bioconductor

# development version:
# remotes::install_github("USERNAME/ARGflanker")
```

## Quick start

```r
library(ARGflanker)

gff <- system.file("extdata", "example.gff3", package = "ARGflanker")
gr  <- classifyGenes(importAnnotations(gff))

# 1. Characterize what flanks each ARG
nb <- characterizeNeighbors(gr, window = 10000)
summarizeNeighbors(nb)

# 2. Learn the dataset-wide signature
neighborhoodProfile(gr)

# 3. Full mobility screen (with optional read-pair bridging)
pairs <- read.csv(
    system.file("extdata", "example_readpairs.csv", package = "ARGflanker"))
res <- MobilityScreen(gr, pairs = pairs, window = 10000)
res
res@scores
```

## Input

- **Gene predictions**: a GFF3 of ORF calls (Prodigal/prokka/bakta style) or a
  `data.frame` with `contig, start, end, strand, gene_id, product`.
- **Annotations**: any product/annotation text; `classifyGenes()` maps it to
  categories. The bundled keyword patterns (`defaultPatterns()`) are a
  placeholder — swap in hits against curated references (CARD/Resfams for ARGs,
  ISfinder/mobileOG for MGE markers) for real work.
- **Read pairs (optional)**: a `data.frame` of paired-end mate mappings
  (`gene_a, category_a, gene_b, category_b`) for cross-contig linkage.

## Status

Early draft (v0.99.0), structured for Bioconductor submission. Function names,
the scoring model and the classification back end are expected to change. The
keyword classifier is illustrative, not validated.

## Development notes

Documentation is written as roxygen2 comments. Regenerate `man/` and `NAMESPACE`
with:

```r
roxygen2::roxygenise()   # or devtools::document()
devtools::check()
```

## License

Artistic-2.0
