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
  df <- data.frame(strings=system(paste("samtools idxstats ", BamFilePath), intern=T))
  
  TotalReads <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% 
    summarise("Total Reads"=sum(as.integer(counts))+sum(as.integer(unaligned)))
  
  TotalAligned <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% 
    summarise("Total Mapped Read Count"=sum(as.integer(counts)))
  
  TotalUnaligned <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% 
    summarise("Total Unmapped Read Count"=sum(as.integer(unaligned)))
  
  out <- data.frame(TotalReads, TotalAligned, TotalUnaligned, "Percent.Mapped" = TotalAligned [[1]] / TotalReads[[1]])
  rownames(out) <- basename(BamFilePath)
  return(out)
}

