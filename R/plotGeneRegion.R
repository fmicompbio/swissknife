#' Prepare GTF file for use with plotGeneRegion
#' 
#' This function sets the names of the transcript and gene ID columns of the
#' gtf file to "transcript" and "gene", removes version tags of the 
#' transcripts/genes and retains only the "exon" entries. The purpose is to 
#' make the file amenable to plotting with Gviz, using the 
#' \code{plotGeneRegion} function.
#' 
#' @param gtf Character scalar, path to gtf file (tested with Ensembl/Gencode
#'   files).
#' @param transcriptIdColumn Character scalar, the column in the gtf file that 
#'   contains the transcript ID.
#' @param geneIdColumn Character scalar, the column in the gtf file that 
#'   contains the gene ID.
#' @param geneSymbolColumn Character scalar, the column in the gtf file that 
#'   contains the gene symbol (if available). Set to \code{""} if not 
#'   available (in which case the gene IDs will be used in its place).
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics subset
#' 
#' @examples
#' gtf <- prepareGTF(gtf = system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
#'                                     package = "swissknife"))
#'                                     
prepareGTF <- function(gtf, transcriptIdColumn = "transcript_id", 
                       geneIdColumn = "gene_id", geneSymbolColumn = "gene_name") {
    .assertPackagesAvailable("rtracklayer")
    gtf <- rtracklayer::import(gtf)
    
    if (!(transcriptIdColumn %in% colnames(S4Vectors::mcols(gtf)))) {
        stop("The gtf file does not have a '", transcriptIdColumn, "' tag")
    }
    if (!(geneIdColumn %in% colnames(S4Vectors::mcols(gtf)))) {
        stop("The gtf file does not have a '", geneIdColumn, "' tag")
    }
    
    ## Set appropriate column names
    idx <- match(c(transcriptIdColumn, geneIdColumn),
                 colnames(S4Vectors::mcols(gtf)))
    colnames(S4Vectors::mcols(gtf))[idx] <- c("transcript", "gene")
    
    ## Strip version numbers from gene and transcript IDs if they exist
    gtf$transcript <- gsub("\\.[0-9]+$", "", gtf$transcript)
    gtf$gene <- gsub("\\.[0-9]+$", "", gtf$gene)
    
    ## Set gene symbol column
    if (!(geneSymbolColumn %in% colnames(S4Vectors::mcols(gtf)))) {
        message("Using the 'gene' column as gene symbols")
        gtf$gene_name <- S4Vectors::mcols(gtf)[["gene"]]
    } else {
        gtf$gene_name <- S4Vectors::mcols(gtf)[[geneSymbolColumn]]
    }
    
    ## Keep only exons
    gtf <- BiocGenerics::subset(gtf, type %in% c("exon", "gene"))
    
    gtf
}

#' Plot gene region
#'
#' Visualize the gene model for a gene of interest, or for all genes in a
#' provided region, and/or show one or more coverage tracks based on bigwig
#' file(s). 
#' 
#' The gene annotation can be provided either as a path to a gtf file, or as a 
#' GRanges object (generated using the \code{prepareGTF} function to ensure 
#' compatibility). The region to display can be determined either by 
#' specifying a gene (ID or symbol) or by specifying a viewing range 
#' (chromosome, start and end positions).
#'
#' @param gtf Character scalar, path to gtf file (tested with Ensembl/Gencode
#'   files).
#' @param granges GRanges object, typically generated from a GTF file using the
#'   \code{prepareGTF} function. This is an alternative to providing the link 
#'   to the gtf file directly, and will take precedence over the \code{gtf} 
#'   argument if provided. 
#' @param chr Character scalar, name of the chromosome to show.
#' @param start,end Numeric scalars, start and end position of the region to
#'   show.
#' @param showgene Character scalar, the gene ID/name to display. Will take 
#'   precedence over positional range specification if provided. 
#' @param bigwigFiles Named character vector, paths to bigwig files.
#' @param bigwigCond Named character vector, the grouping of the bigwig files
#'   (used for coloring of the coverage tracks).
#' @param geneTrackTitle Character scalar, name of the gene track.
#' @param transcriptIdColumn Character scalar, the column in the gtf file that 
#'   contains the transcript ID. Passed to \code{prepareGTF}.
#' @param geneIdColumn Character scalar, the column in the gtf file that 
#'   contains the gene ID. Passed to \code{prepareGTF}.
#' @param geneSymbolColumn Character scalar, the column in the gtf file that 
#'   contains the gene symbol (if available). Set to \code{""} if not 
#'   available (in which case the gene IDs will be used in its place).
#'   Passed to \code{prepareGTF}.
#' @param lowerPadding,upperPadding Numeric scalars, setting the amount of 
#'   padding in the lower and upper range of the plot, respectively. For 
#'   example, a value of 0.05 will expand the range by 
#'   0.05 * (max coordinate - min coordinate) in the specified direction. 
#' @param colorByStrand Logical scalar, determining whether gene features are 
#'   colored by the annotated strand.
#' @param featureColors Named character vector of length 4, with elements 
#'   \code{plusmain}, \code{minusmain}, \code{plusother}, \code{minusother}, 
#'   giving the colors to use for the features if \code{colorByStrand} is TRUE. 
#' @param condColors Either NULL or a named character vector (with the same 
#'   names as the unique values of \code{bigwigCond}), giving the colors to use
#'   for the coverage tracks if \code{bigwigCond} is provided. 
#' @param scaleDataTracks Logical scalar, indicating whether the data tracks 
#'   should be scaled to have the same y-axis limits. 
#' @param plotTitle Character scalar, the title of the final plot. If 
#'   \code{NULL} (the default), it will be automatically defined based on 
#'   the displayed gene or region.
#' @param ... Additional arguments to be passed to \code{Gviz::plotTracks}.
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom BiocGenerics subset start end strand
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges overlapsAny IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors %in%
#' @importFrom methods is
#'   
#' @examples 
#' if (requireNamespace("Gviz", quietly = TRUE)) {
#'     gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf", 
#'                            package = "swissknife")
#'     plotGeneRegion(gtf = gtffile, 
#'                    showgene = "Tnfaip3")
#'                
#'     bwf <- system.file("extdata/plotGeneRegion/mnase_mm10.bw", 
#'                        package = "swissknife")
#'     names(bwf) <- "bwf1"
#'     plotGeneRegion(gtf = gtffile, 
#'                    bigwigFiles = bwf,
#'                    chr = "chr10", start = 20000000, end = 20005000)
#'     plotGeneRegion(bigwigFiles = bwf,
#'                    chr = "chr10", start = 20000000, end = 20005000)
#'                
#'     bwf2 <- c(bwf, bwf)
#'     names(bwf2) <- c("bwf1", "bwf2")
#'     bwc2 <- c("c1", "c2")
#'     names(bwc2) <- names(bwf2)
#'     plotGeneRegion(gtf = gtffile, bigwigFiles = bwf2, bigwigCond = bwc2, 
#'                    showgene = "Map3k5")
#' }
#' 
plotGeneRegion <- function(gtf = "", granges = NULL, chr = "", 
                           start = NA_real_, end = NA_real_, showgene = "", 
                           bigwigFiles = "", bigwigCond = "",
                           geneTrackTitle = "Genes", 
                           transcriptIdColumn = "transcript_id", 
                           geneIdColumn = "gene_id",
                           geneSymbolColumn = "gene_name",
                           lowerPadding = 0.15, 
                           upperPadding = 0.05,
                           colorByStrand = FALSE, 
                           featureColors = c(plusmain = "#0E14D0", 
                                             minusmain = "#D0350E",
                                             plusother = "#9E9BEB", 
                                             minusother = "#DA907E"),
                           condColors = NULL,
                           scaleDataTracks = FALSE, 
                           plotTitle = NULL, ...) {
    options(ucscChromosomeNames = FALSE)
    
    ## ---------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## ---------------------------------------------------------------------- ##
    .assertPackagesAvailable(c("Gviz", "rtracklayer"))
    ## Check that input classes are correct
    if (!methods::is(gtf, "character") || length(gtf) != 1) {
        stop("'gtf' must be a character scalar")
    }
    if (!is.null(granges) && !methods::is(granges, "GRanges")) {
        stop("'granges' must be either NULL or a GRanges object")
    }
    if (!methods::is(chr, "character") || length(chr) != 1) {
        stop("'chr' must be a character scalar")
    }
    if (length(start) != 1 || length(end) != 1) {
        stop("'start' and 'end' must be numeric scalars")
    }
    if ((!is.na(start) && !methods::is(start, "numeric")) || 
        (!is.na(end) && !methods::is(end, "numeric"))) {
        stop("'start' and 'end' must be numeric scalars")
    }
    if (!is.na(start) && !is.na(end) && start >= end) {
        stop("'end' must be strictly larger than 'start'")
    }
    if (!methods::is(showgene, "character") || length(showgene) != 1) {
        stop("'showgene' must be a character scalar")
    }
    if (!methods::is(bigwigFiles, "character")) {
        stop("'bigwigFiles' must be a character vector")
    }
    if (any(bigwigFiles != "") && any(is.null(names(bigwigFiles)))) {
        stop("'bigwigFiles' must be a named character vector")
    }
    if (!methods::is(bigwigCond, "character")) {
        stop("'bigwigCond' must be a character vector")
    }
    if (any(bigwigCond != "") && any(is.null(names(bigwigCond)))) {
        stop("'bigwigCond' must be a named character vector")
    }
    if (any(bigwigFiles != "") && any(bigwigCond != "") && 
        any(names(bigwigFiles) != names(bigwigCond))) {
        stop("'bigwigFiles' and 'bigwigCond' must have the same names")
    }
    if (!methods::is(geneTrackTitle, "character") || 
        length(geneTrackTitle) != 1) {
        stop("'geneTrackTitle' must be a character scalar")
    }
    if (!methods::is(lowerPadding, "numeric") || 
        !methods::is(upperPadding, "numeric") || 
        length(lowerPadding) != 1 || length(upperPadding) != 1) {
        stop("'lowerPadding' and 'upperPadding' must be numeric scalars")
    }
    if (!is.logical(colorByStrand) || 
        length(colorByStrand) != 1) {
        stop("'colorByStrand' must be a logical scalar")
    }
    if (!methods::is(featureColors, "character") || 
        length(featureColors) != 4) {
        stop("'featureColors' must be a character vector of length 4")
    }
    if (length(intersect(names(featureColors), 
                         c("plusmain", "minusmain", "plusother", "minusother"))) != 4) {
        stop("'featureColors must have elements 'plusmain', 'minusmain', ",
             "'plusother', 'minusother'")
    }
    if (!is.null(condColors)) {
        if (!is(condColors, "character")) {
            stop("'condColors' must be a character vector or NULL")
        }
        if (any(bigwigCond != "") && 
            length(intersect(unique(bigwigCond), names(condColors))) != 
            length(condColors)) {
            stop("'condColors' must be a named vector with names equal to the ",
                 "values of bigwigCond")
        }
    }
    if (length(scaleDataTracks) != 1 || !is.logical(scaleDataTracks)) {
        stop("'scaleDataTracks' must be a logical scalar")
    }
    if (!is.null(plotTitle)) {
        if (length(plotTitle) != 1 || !is.character(plotTitle)) {
            stop("'plotTitle' must be a character scalar")
        }
    }
    ## Must have at least one of bigwigFiles, gtf and granges
    if (all(bigwigFiles == "") && is.null(granges) && gtf == "") {
        stop("'At least one of 'bigwigFiles', 'granges' and 'gtf' must be provided")
    }
    
    ## Read gtf file if provided
    if (is.null(granges) && gtf != "") {
        granges <- prepareGTF(gtf, transcriptIdColumn = transcriptIdColumn,
                              geneIdColumn = geneIdColumn,
                              geneSymbolColumn = geneSymbolColumn)
    }
    
    ## If granges file is not given, the viewing region must be set
    if (is.null(granges) && (chr == "" || is.na(start) || is.na(end))) {
        stop("Either an annotation has to be given, or the viewing region must be set")
    }
    
    ## Extract full gene loci into its own GRanges object, in order to 
    ## show also features where only the introns overlap the selected region
    if (!is.null(granges)) {
        grangesgene <- BiocGenerics::subset(granges, type == "gene")
        granges <- BiocGenerics::subset(granges, type == "exon")
    }
    
    if (!is.na(start)) {
        start <- as.integer(start)
    }
    if (!is.na(end)) {
        end <- as.integer(end)
    }
    
    ## ---------------------------------------------------------------------- ##
    ## Prepare the annotation
    ## ---------------------------------------------------------------------- ##
    if (!is.null(granges)) {
        ## Strip version number from the gene of interest if it exists
        showgene <- gsub("\\.[0-9]+$", "", showgene)
        
        if (showgene == "" && (chr == "" || is.na(start) || is.na(end))) {
            stop("Either a gene has to be given, or the viewing region must be set")
        }
        
        ## If a gene has been defined, set the viewing range accordingly
        if (showgene != "") {
            gr <- BiocGenerics::subset(granges, tolower(gene) == tolower(showgene) | 
                                           tolower(gene_name) == tolower(showgene))
            ## Select only one gene if there are many with the same name
            if (length(unique(gr$gene)) > 1) {
                warning("Multiple gene IDs with the same name, choosing the first")
            }
            gr <- BiocGenerics::subset(gr, gene == gene[1])
            chr <- unique(GenomeInfoDb::seqnames(gr))
            start <- min(BiocGenerics::start(gr))
            end <- max(BiocGenerics::end(gr))
        } else {
            gr <- grangesgene[IRanges::overlapsAny(
                grangesgene,
                GenomicRanges::GRanges(seqnames = chr,
                                       ranges = IRanges::IRanges(start = start,
                                                                 end = end),
                                       strand = "*")), ]
            gr <- BiocGenerics::subset(granges, gene %in% gr$gene)
        }
        
        ## Other features in the region
        grogene <- BiocGenerics::subset(grangesgene, !(gene %in% gr$gene))
        gro <- grogene[IRanges::overlapsAny(
            grogene,
            GenomicRanges::GRanges(seqnames = chr,
                                   ranges = IRanges::IRanges(start = start,
                                                             end = end),
                                   strand = "*"))]
        gro <- BiocGenerics::subset(granges, gene %in% gro$gene)
        
        grtr <- Gviz::GeneRegionTrack(gr, showId = TRUE, col = NULL, fill = "gray80",
                                      name = geneTrackTitle, col.title = "black")
        
        if (colorByStrand) {
            grtr2 <- Gviz::GeneRegionTrack(gro, showId = TRUE, col = NULL, 
                                           fill = "white",
                                           name = "", col.title = "black")
            if (length(grtr) > 0) {
                Gviz::feature(grtr) <- ifelse(BiocGenerics::strand(grtr) == "+", 
                                              "plusmain", "minusmain")
            }
            if (length(grtr2) > 0) {
                Gviz::feature(grtr2) <- ifelse(BiocGenerics::strand(grtr2) == "+", 
                                               "plusother", "minusother")
            }
        } else {
            grtr2 <- Gviz::GeneRegionTrack(gro, showId = TRUE, col = "black", 
                                           fill = "white",
                                           name = "", col.title = "black")
        }
    } else {
        gr <- gro <- grtr <- grtr2 <- NULL
    }    
    
    ## ---------------------------------------------------------------------- ##
    ## Set title (if not provided) and viewing range
    ## ---------------------------------------------------------------------- ##
    if (is.null(plotTitle)) {
        ## Define the title for the plot
        if (showgene != "" && !is.null(gr)) {
            if (all(gr$gene == gr$gene_name)) {
                plotTitle <- unique(gr$gene)
            } else {
                plotTitle <- unique(paste0(gr$gene, " (", gr$gene_name, ")"))
            }
        } else {
            plotTitle <- paste0(chr, ":", start, "-", end)
        }
    }
    
    ## Set min and max coord for the plot (add some padding to each side)
    minCoord <- start - lowerPadding*(end - start)
    maxCoord <- end + upperPadding*(end - start)
    
    ## ---------------------------------------------------------------------- ##
    ## Prepare bigWig files
    ## ---------------------------------------------------------------------- ##
    ## Reformat bigWig file paths and names (provided to the function as 
    ## character strings)
    if (!any(bigwigFiles == "")) {
        
        ## Determine the y limits
        ## Thanks to http://adomingues.github.io/2016/11/13/max-coverage-in-bigwigs/
        ## for the idea
        if (!scaleDataTracks) {
            ylim <- NULL
        } else {
            g <- GenomicRanges::GRanges(seqnames = chr,
                                        ranges = IRanges(start = minCoord,
                                                         end = maxCoord),
                                        strand = "*")
            ymax <- max(vapply(bigwigFiles, function(f) {
                ff <- rtracklayer::import(f, format = "bw", which = g)
                if (length(ff) > 0) {
                    return(max(ff$score))
                } else {
                    return(0)
                }
            }, 0))
            ylim <- c(0, ymax)
        }
        
        ## ---------------------------------------------------------------------- ##
        ## Define colors if bigwigCond is provided
        ## ---------------------------------------------------------------------- ##
        if (!is.null(condColors)) {
            usecol <- condColors[match(bigwigCond, names(condColors))]
        } else {
            ## Define colors for coverage tracks
            color_list <- rep(c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "#F7EE55",
                                "#90C987", "#777777", "#E8601C", "#1965B0", "#882E72",
                                "#F6C141", "#4EB265", "#CAEDAB"), 
                              ceiling(length(unique(bigwigCond))/13))
            
            if (length(bigwigCond) > 1 || bigwigCond != "") {
                usecol <- color_list[match(bigwigCond, 
                                           unique(bigwigCond))]
            } else {
                usecol <- rep("gray", length(bigwigFiles))
            }
            names(usecol) <- names(bigwigCond)
        }
        
        ## ------------------------------------------------------------------ ##
        ## Prepare final plot
        ## ------------------------------------------------------------------ ##
        ## Set up coverage tracks
        tracks <- lapply(seq_along(bigwigFiles), function(i) {
            assign(paste0("covtr", i), 
                   Gviz::DataTrack(range = bigwigFiles[i],
                                   type = "histogram",
                                   name = names(bigwigFiles)[i],
                                   col.title = "black",
                                   fill = usecol[i],
                                   col = usecol[i],
                                   col.histogram = usecol[i],
                                   fill.histogram = usecol[i],
                                   ylim = ylim))
        })
    } else {
        tracks <- NULL
    }
    
    ## Add genome axis track
    if (length(grtr) == 0) grtr <- NULL
    if (length(grtr2) == 0) grtr2 <- NULL
    tracks <- c(tracks, Gviz::GenomeAxisTrack(), grtr, grtr2)
    
    ## Plot tracks
    Gviz::plotTracks(tracks, chromosome = chr, from = minCoord, 
                     to = maxCoord, main = plotTitle, 
                     transcriptAnnotation = "transcript",
                     min.width = 0, min.distance = 0, collapse = FALSE,
                     plusmain = featureColors["plusmain"], 
                     minusmain = featureColors["minusmain"],
                     plusother = featureColors["plusother"], 
                     minusother = featureColors["minusother"], ...)
}
