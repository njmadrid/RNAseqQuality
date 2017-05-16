#' Get a summary of the unique reads from bam file. Duplicate reads are defined as having identical pos and cigar values.
#' 
#' unique aligned reads = total number of reads that are not secondary;
#' uniquely aligned reads without duplicates = total number of reads that are not secondary and not duplicates;
#' duplicate reads = total number of reads that are duplicates;
#'
#' @param BamFilePath location of aligned bam file. Corresponding bai index file must be in the same directory.
#' @keywords SummarizeUniqueReads
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame containing a summary of the unique reads from bam file. Duplicate reads are defined as having identical pos and cigar values.
#' @examples
#' BamFilePaths <- dir(path="/home/ubuntu/Projects/RNAseqPackage/bams", pattern=".*.bam$", full.names=T)
#' UniqueReadsDF = SummarizeUniqueReads(BamFilePaths[1])
#' @export
SummarizeUniqueReads <- function(BamFilePath) {
  # unique aligned reads = total number of reads that are not secondary
  # uniquely aligned reads without duplicates = total number of reads that are not secondary and not duplicates
  # duplicate reads = total number of reads that are duplicates
  
  param <- ScanBamParam(what=c("flag", "rname", "pos", "cigar"), flag=scanBamFlag(isUnmappedQuery=F))
  ReadsDF <- as.data.frame(scanBam(BamFilePath, param=param)[[1]])
  ReadsDF$isDuplicate <- duplicated(ReadsDF[ , -which(colnames(ReadsDF) == "flag")])
  # flag = 256 designates a secondary alignment
  UniqueAlignedReads <- nrow(ReadsDF[ReadsDF$flag != 256, ])
  UniquelyAlignedReadsWithoutDuplicates <- nrow(ReadsDF[ReadsDF$flag != 256 & ReadsDF$isDuplicate == F, ])
  DuplicateReads <- nrow(ReadsDF[which(ReadsDF$isDuplicate == T), ])
  
  out <- data.frame("unique aligned reads" = UniqueAlignedReads,
                    "uniquely aligned reads without duplicates" = UniquelyAlignedReadsWithoutDuplicates,
                    "duplicate reads" = DuplicateReads)
  
  return(out)
}