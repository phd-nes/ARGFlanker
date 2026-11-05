# ARGflanker 0.99.0

* Initial draft for Bioconductor submission.
* Core neighbourhood characterization: `characterizeNeighbors()`,
  `summarizeNeighbors()`, `neighborhoodProfile()`.
* Mobility evidence and scoring: `screenMobility()`, `linkByReadPairs()`,
  `mobilityScore()`, and the `MobilityScreen` class / constructor.
* I/O and classification: `importAnnotations()`, `classifyGenes()`,
  `defaultPatterns()`.
* Visualisation: `plotNeighborhood()`.
* Bundled example data (`inst/extdata/example.gff3`,
  `inst/extdata/example_readpairs.csv`) and a vignette.
