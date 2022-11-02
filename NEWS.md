# swissknife 0.39

* Add function to set LD_LIBRARY_PATH in conda environments

# swissknife 0.38

* Add plotTitle argument to plotGeneRegion

# swissknife 0.37

* Add ellipsis to plotGeneRegion to allow additional arguments to plotTracks

# swissknife 0.36

* Add .assertPackagesAvailable utility function and unit test to templates
* Update internally used .assertPackagesAvailable to be identical to the ones in templates

# swissknife 0.35

* Expand capabilities of utility functions, add unit tests
* fix bowtie executable location logic for windows

# swissknife 0.34

* Bug fix - expand documentation in utility functions.

# swissknife 0.33

* Add addUtilsFunctions to include utility functions in other packages

# swissknife 0.32

* Add getInsertSizeDistFromBam to tabulate insert sizes from paired-end alignments in bam files

# swissknife 0.31

* Fix readSampleTsvs so that provided patterns must match from the beginning of the filename

# swissknife 0.30

* Speed-up package load time by moving rarely used Depends to Suggests

# swissknife 0.29

* Refactor getGenomicTiles, adding a new annotateRegions that accepts any GRanges as input

# swissknife 0.28

* Truncate chromosome names after first whitespace in getGenomicTiles when checking consistency

# swissknife 0.27

* Add getGenomicTiles to obtain annotated regions tiling a given genome

# swissknife 0.26

* speed-up getMappableRegions by dumping k-mers in C++

# swissknife 0.25

* fix getMappableRegions to work for fully mappable chromosomes

# swissknife 0.24

* Add getMappableRegions to obtain mappable regions in a genome

# swissknife 0.23

* Add specificityScore to calculate (tissue) specificity of gene expression

# swissknife 0.22

* ignore NA values in normGenesetExpression

# swissknife 0.21

* Add labelCells to label cells using sets of marker genes
* Add normGenesetExpression to calculate normalized expression of a set of genes

# swissknife 0.20

* set seed in unit test for sampleControlElements

# swissknife 0.19

* fix coloring bugin in plotGeneRegion

# swissknife 0.18

* setup github actions and remove travis

# swissknife 0.17

* Add R/selVarGenes.R to select variable genes in a single cell RNA-seq experiment
* configured travis/github to automatically run pkgdown::deploy_site_github()

# swissknife 0.16

* Show features that overlap selected regions only with their introns in plotGeneRegion

# swissknife 0.15

* Allow specification of colors for coverage tracks in plotGeneRegion

# swissknife 0.14

* Allow coloring features by strand in plotGeneRegion

# swissknife 0.13

* Add parsePkgVersions
* Add plotGeneRegion

# swissknife 0.12

* Add col2hex

# swissknife 0.10 / 0.11

* Add ypixels, pointsize argument to plotBitScatter

# swissknife 0.9

* Add R/valueToColor.R

# swissknife 0.8

* plotBitScatter(x, y = NULL) now also works for x being a two-column matrix

# swissknife 0.7

* add functions to calculate and plot phasograms, e.g. to measure nucleosome
  repeat length in MNase data (R/phasograms.R, src/phasograms.cpp)

# swissknife 0.6

* added tests for sampleControlElements() to cover all code (coverage gaps detected
  with covr::zero_coverage(covr::package_coverage()))

# swissknife 0.5

* sampleControlElements() now catches cases with too few control elements to match the target distribution

# swissknife 0.4

* sampleControlElements(x = list(), ...) now works correctly

# swissknife 0.3

* Add R/readSampleTsvs.R:readSampleTsvs()

# swissknife 0.2

* Add R/plotting.R:plotBitScatter()

# swissknife 0.1

* Initial version
* add R/misc.R:sampleControlElements()
