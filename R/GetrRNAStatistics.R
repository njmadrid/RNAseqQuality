#' Get the total number of rRNA read counts from bam file.
#' samtools idxstats must be installed on the local system.
#' 
#' Note that samtools idxstats may or may not include secondary alignments. This will depend on how the reads
#' were aligned. See the documention for samtools idxstats and the aligner used to create the bam file.
#' The bam file must have been aligned to a reference genome containing chromosone GL000220.1
#' total reads mapping to GL0000220.1 or know rRNA regions defined by Ensembl
#' 
#' @param BamFilePath location of aligned bam file. Corresponding bai index file must be in the same directory.
#' @keywords GetrRNAStatistics
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return An integer for the total number of rRNA read counts in a given bam file.
#' @examples
#' BamFilePaths <- dir(path="/home/ubuntu/Projects/RNAseqPackage/bams", pattern=".*.bam$", full.names=T)
#' # samtools idxstats must be installed on the local system.
#' # Note that samtools idxstats may or may not include secondary alignments. This will depend on how the reads
#' # were aligned. See the documention for samtools idxstats and the aligner used to create the bam file.
#' rRNAStatsDF = GetrRNAStatistics(BamFilePaths[1])
#' @export
GetrRNAStatistics <- function(BamFilePath) {
  df <- data.frame(strings=system(paste("samtools idxstats ", BamFilePath), intern = T))
  
  TotalAligned <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% 
    summarise("Total Alignment Count"=sum(as.integer(counts)))
  
  TotalUnaligned <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% 
    summarise("Total Unmapped Read Count"=sum(as.integer(unaligned)))
  
  rrna <- df %>% separate(strings, into=c("gene","length","counts","unaligned"),sep = "\t") %>% filter(gene=="GL000220.1") %>% select(counts)
  out <- data.frame(TotalAligned, TotalUnaligned, "Alignments Mapped to rRNA"=rrna[[1]])
  rownames(out) <- basename(BamFilePath)
  return(out)
}