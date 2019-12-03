#' Prepare GTF file for use with plotGeneRegion
#' 
#' @param gtf Character scalar, path to gtf file (tested with Ensembl/Gencode files)
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics subset
#' 
prepareGTF <- function(gtf) {
    gtf <- rtracklayer::import(gtf)
    
    ## Set appropriate column names
    idx <- match(c("transcript_id", "gene_id"),
                 colnames(S4Vectors::mcols(gtf)))
    colnames(S4Vectors::mcols(gtf))[idx] <- c("transcript", "gene")
    if (!("gene_name" %in% colnames(S4Vectors::mcols(gtf)))) {
        gtf$gene_name <- gtf$gene
    }
    
    ## Keep only exons
    gtf <- BiocGenerics::subset(gtf, type == "exon")
    
    ## Strip version numbers from gene and transcript IDs if they exist
    gtf$transcript <- gsub("\\.[0-9]+$", "", gtf$transcript)
    gtf$gene <- gsub("\\.[0-9]+$", "", gtf$gene)
    
    gtf
}

#' Plot gene region
#'
#' Visualize the gene model for a gene of interest, or for all genes in a
#' provided region. Also show one or more coverage tracks based on bigwig
#' file(s).
#'
#' @param gtf Character scalar, path to gtf file (tested with Ensembl/Gencode
#'   files)
#' @param granges GRanges object, typically generated from a GTF file using the
#'   \code{prepareGTF} function
#' @param chr Character scalar, name of chromosome to show
#' @param start,end Numeric scalars, start and end position of the region to
#'   show
#' @param showgene Character scalar, the gene ID/name to display
#' @param bigwigFiles Named character vector, paths to bigwig files
#' @param bigwigCond Named character vector, the grouping of the bigwig files
#'   (used for coloring of the coverage tracks)
#' @param geneTrackTitle Character scalar, name of the gene track
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom BiocGenerics subset start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges overlapsAny IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom Gviz GeneRegionTrack DataTrack GenomeAxisTrack plotTracks
#' @importFrom S4Vectors %in%
#' @importFrom methods is
#'   
#' @examples 
#' gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf", 
#'                        package = "swissknife")
#' plotGeneRegion(gtf = gtffile, 
#'                showgene = "Tnfaip3")
#'                
#' bwf <- system.file("extdata/plotGeneRegion/mnase_mm10.bw", package = "swissknife")
#' names(bwf) <- "bwf1"
#' plotGeneRegion(gtf = gtffile, 
#'                bigwigFiles = bwf,
#'                chr = "chr10", start = 20000000, end = 20005000)
#'                
#' bwf2 <- c(bwf, bwf)
#' names(bwf2) <- c("bwf1", "bwf2")
#' bwc2 <- c("c1", "c2")
#' names(bwc2) <- names(bwf2)
#' plotGeneRegion(gtf = gtffile, bigwigFiles = bwf2, bigwigCond = bwc2, 
#'                showgene = "Map3k5")
#' 
plotGeneRegion <- function(gtf = "", granges = NULL, chr = "", 
                           start = NA_real_, end = NA_real_, showgene = "", 
                           bigwigFiles = "", bigwigCond = "",
                           geneTrackTitle = "Genes") {
    options(ucscChromosomeNames = FALSE)
    
    ## ---------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## ---------------------------------------------------------------------- ##
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
    if (!methods::is(start, "numeric") || !methods::is(end, "numeric") || 
        length(start) != 1 || length(end) != 1) {
        stop("'start' and 'end' must be numeric scalars")
    }
    if (!methods::is(showgene, "character") || length(showgene) != 1) {
        stop("'showgene' must be a character scalar")
    }
    if (!methods::is(bigwigFiles, "character")) {
        stop("'bigwigFiles' must be a character vector")
    }
    if (bigwigFiles != "" && any(is.null(names(bigwigFiles)))) {
        stop("'bigwigFiles' must be a named character vector")
    }
    if (!methods::is(bigwigCond, "character")) {
        stop("'bigwigCond' must be a character vector")
    }
    if (bigwigCond != "" && any(is.null(names(bigwigCond)))) {
        stop("'bigwigCond' must be a named character vector")
    }
    if (bigwigFiles != "" && bigwigCond != "" && 
        names(bigwigFiles) != names(bigwigCond)) {
        stop("'bigwigFiles' and 'bigwigCond' must have the same names")
    }
    if (!methods::is(geneTrackTitle, "character") || 
        length(geneTrackTitle) != 1) {
        stop("'geneTrackTitle' must be a character scalar")
    }
    
    ## Must have at least one of bigwigFiles, gtf and granges
    if (bigwigFiles == "" && is.null(granges) && gtf == "") {
        stop("'At least one of 'bigwigFiles', 'granges' and 'gtf' must be provided")
    }
    
    ## Read gtf file if provided
    if (is.null(granges) && gtf != "") {
        granges <- prepareGTF(gtf)
    }
    
    ## If granges file is not given, the viewing region must be set
    if (is.null(granges) && (chr == "" || is.na(start) || is.na(end))) {
        stop("Either an annotation has to be given, or the viewing region must be set")
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
            gr <- granges[IRanges::overlapsAny(
                granges,
                GenomicRanges::GRanges(seqnames = chr,
                                       ranges = IRanges::IRanges(start = start,
                                                                 end = end),
                                       strand = "*")), ]
        }
        
        ## Other features in the region
        gro <- granges[IRanges::overlapsAny(
            granges,
            GenomicRanges::GRanges(seqnames = chr,
                                   ranges = IRanges::IRanges(start = start,
                                                             end = end),
                                   strand = "*"))]
        gro <- gro[!(S4Vectors::`%in%`(gro, gr))]
        
        grtr <- Gviz::GeneRegionTrack(gr, showId = TRUE, col = NULL, fill = "gray80",
                                      name = geneTrackTitle, col.title = "black")
        grtr2 <- Gviz::GeneRegionTrack(gro, showId = TRUE, col = "black", fill = "white",
                                       name = "", col.title = "black")
    } else {
        gr <- gro <- grtr <- grtr2 <- NULL
    }    
    
    ## ---------------------------------------------------------------------- ##
    ## Set title and viewing range
    ## ---------------------------------------------------------------------- ##
    ## Define the title for the plot
    if (showgene != "" && !is.null(gr)) {
        if (all(gr$gene == gr$gene_name)) {
            plot_title <- unique(gr$gene)
        } else {
            plot_title <- unique(paste0(gr$gene, " (", gr$gene_name, ")"))
        }
    } else {
        plot_title <- paste0(chr, ":", start, "-", end)
    }
    
    ## Set min and max coord for the plot (add some padding to each side)
    minCoord <- start - 0.15*(end - start)
    maxCoord <- end + 0.05*(end - start)
    
    ## ---------------------------------------------------------------------- ##
    ## Prepare bigWig files
    ## ---------------------------------------------------------------------- ##
    ## Reformat bigWig file paths and names (provided to the function as 
    ## character strings)
    if (bigwigFiles != "") {
        
        ## ---------------------------------------------------------------------- ##
        ## Define colors if bigwigCond is provided
        ## ---------------------------------------------------------------------- ##
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
                                   fill.histogram = usecol[i]))
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
                     to = maxCoord, main = plot_title, 
                     transcriptAnnotation = "transcript",
                     min.width = 0, min.distance = 0, collapse = FALSE)
}
