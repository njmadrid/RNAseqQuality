#' Overlap bedgraph coverage file (*.bedgraph.gz) with known editing sites (e.g. 
#' from RADAR: Rigorously Annotated Database of A-to-I RNA editing or from DARNED: 
#' a DAtabase of RNa EDiting in humans) to obtain the number of RNA editing sites with coverage.
#' 
#' @param BedGraphFilePath location of bedgraph coverage file (*.bedgraph.gz) created via bedtools genomecov -bga -split.
#' @param EditsGRanges GRanges object of known RNA edit sites.
#' @keywords OverlapEditSites
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame containing the mcols data from a Hits object returned in findOverlaps.
#' @examples
#' radarFile <- "dependencies/Human_AG_all_hg19_v2.txt"
#' radar <- read.table(radarFile, header=T)
#' 
#' # fix contig names
#' radar$chromosome <- sub("chr","",radar$chromosome)
#' radar$chromosome[radar$chromosome=="M"] <- "MT" 
#' radarRanges <- with(radar, GRanges(seq=chromosome, ranges=IRanges(start=position,end=position), strand=strand, gene=gene,
#'                                    annot1=annot1, annot2=annot2, alu=alu., non_alu_repetitive=non_alu_repetitive.,
#'                                    conservation_chimp=conservation_chimp, conservation_rhesus,conservation_rhesus,
#'                                    conservation_mouse=conservation_mouse))
#' 
#' # sort so the seqlevels actually match
#' radarRanges <- sortSeqlevels(radarRanges)
#' radarRanges <- sort(radarRanges)
#' # Get file paths of bedgraph files
#' BedGraphFiles <- dir(path="/home/ubuntu/Projects/RNAseqPackage/bedgraphs", full.names=T)
#' RadarOverlapsDF <- do.call("cbind", mclapply(BedGraphFiles, function(f) OverlapEditSites(f, radarRanges), mc.cores=11))
#' colnames(RadarOverlapsDF) <- basename(BedGraphFiles)
#' RadarOverlapsDF <- cbind(as.data.frame(mcols(radarRanges)), as.data.frame(RadarOverlapsDF))
#' mcols(radarRanges) <- RadarOverlapsDF
#' @export
OverlapEditSites <- function(BedGraphFilePath, EditsGRanges) {
  bedgraph <- import.bedGraph(BedGraphFilePath)
  scores <- mcols(bedgraph[findOverlaps(EditsGRanges, bedgraph, select='first'),])$score
  return(scores)
}

