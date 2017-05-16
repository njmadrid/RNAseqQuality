#' Get a summary of alignment statistics for an aligned bam file as defined by samtools idxstats.
#'
#' @param BamFilePath location of aligned bam file. Corresponding bai index file must be in the same directory.
#' @keywords GetAlignmentStatistics
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame of a summary of alignment statistics for an aligned bam file as defined by samtools idxstats.
#' @examples
#' BamFilePaths <- dir(path="/home/ubuntu/Projects/RNAseqPackage/bams", pattern=".*.bam$", full.names=T)
#' AlignmentStatsDF = GetAlignmentStatistics(BamFilePaths[1])
#' @export
GetAlignmentStatistics <- function(BamFilePath) {
  TextDF <- data.frame(strings=system(paste("samtools idxstats ", BamFilePath), intern=T))
  TotalReads <- tidyr::separate(TextDF, strings, into=c("gene", "length", "counts", "unaligned"), sep = "\t")
  TotalReads <- plyr::summarise(TotalReads, "Total Reads"=sum(as.integer(counts))+sum(as.integer(unaligned)))
  
  TextDF <- data.frame(strings=system(paste("samtools idxstats ", BamFilePath), intern=T))
  TotalAligned <- tidyr::separate(TextDF, strings, into=c("gene", "length", "counts", "unaligned"), sep = "\t")
  TotalAligned <- plyr::summarise(TotalAligned, "Total Mapped Read Count"=sum(as.integer(counts)))
  
  TextDF <- data.frame(strings=system(paste("samtools idxstats ", BamFilePath), intern=T))
  TotalUnaligned <- tidyr::separate(TextDF, strings, into=c("gene", "length", "counts", "unaligned"), sep = "\t")
  TotalUnaligned <- plyr::summarise(TotalUnaligned, "Total Unmapped Read Count"=sum(as.integer(unaligned)))
  
  out <- data.frame(TotalReads, TotalAligned, TotalUnaligned, "Percent.Mapped" = TotalAligned[[1]] / TotalReads[[1]] * 100)
  rownames(out) <- basename(BamFilePath)
  return(out)
}

