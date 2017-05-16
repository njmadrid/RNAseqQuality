#' Calculate a coverage score based on read depth and also gene length (i.e. number of coding nucleotides).
#' See the Value section of this help page to see how the coverage score is calculated. 
#'
#' @param BedGraphRefSeqOverlapsDF list (1 data frame per file) of overlap scores for a bedgraph coverage file (*.bedgraph.gz) compared to RefSeq defined gene regions.
#' @param TranscriptLengths lengths (i.e. number of coding nucleotides) of each gene defined in RefSeq.
#' @keywords GetTranscriptomeCoverage
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame containing a coverage score for each gene. The coverage score is calculated as follows:
#' @return Coverage Score = sum(width * score) / TranscriptLengths;
#' @return where:
#' @return score = overlap score per exon as determined by GenomicRanges findOverlaps in GetBedGraphRefSeqOverlaps;
#' @return width = width of exon (i.e. number of nucleotides);
#' @return TranscriptLengths = lengths (i.e. number of coding nucleotides) of each gene defined in RefSeq;
#' @examples
#' refseq <- import.gff("RefSeq_hg19_exons_nodups_021114g1k.gff")
#' refseq$gene <- sub("_exon_[0-9]+", "", refseq$gene_id)
#' refseq$width <- width(refseq)
#' txlength <- tapply(refseq$width, refseq$gene, sum)
#' 
#' # Get file paths of aligned bedgraph files
#' BedGraphFiles <- dir(path="/home/ubuntu/Projects/RNAseqPackage/bedgraphs", full.names=T)
#' 
#' # coverage calculated w/ bedtools genomecov -bga -split
#' BedgraphList <- lapply(BedGraphFiles[1:3], function(bg) import.bedGraph(bg))
#' BedGraphRefSeqOverlapsList <- lapply(BedgraphList, function(bg) GetBedGraphRefSeqOverlaps(bg, refseq))
#' BedGraphRefSeqOverlapsDFList <- lapply(BedGraphRefSeqOverlapsList, function(bg) as.data.frame(bg@elementMetadata@listData))
#' for (i in 1:length(BedGraphRefSeqOverlapsDFList)) { BedGraphRefSeqOverlapsDFList[[i]]$width <- BedGraphRefSeqOverlapsList[[i]]@ranges@width }
#' CoverageDFList <- lapply(BedGraphRefSeqOverlapsDFList, function(bg) GetTranscriptomeCoverage(bg, txlength))
#' avgReadDepthPerGene <- unlist(lapply(CoverageDFList, function(df) mean(df$score)))
#' mmdOverTranscriptome <- unlist(lapply(CoverageDFList, function(df) sum(df$TotalCoverage) / sum(df$TranscriptLengths)))
#' @export
GetTranscriptomeCoverage <- function(BedGraphRefSeqOverlapsDF, TranscriptLengths) {
  df <- BedGraphRefSeqOverlapsDF
  df <- ddply(df, .(gene), summarise, TotalCoverage = sum(width * score))
  df$TranscriptLengths <- txlength[as.character(df$gene)]
  df$score <- df$TotalCoverage / unname(df$TranscriptLengths)
  return(df)
}
