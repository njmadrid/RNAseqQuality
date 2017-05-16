#' Read and organize htseq *_summary_stats.txt files into one data structure.
#'
#' @param HtseqFilePaths htseq *_summary_stats.txt file paths.
#' @keywords GetHtseqSummaryStats
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame containing the summary statistics from each HtSeq file in HtseqFilePaths.
#' @examples
#' HtseqFilePaths <- dir(path=".", pattern=".*_summary_stats.txt", full.names=T)
#' HtseqStatsPerFile <- GetHtseqSummaryStats(HtseqFilePaths)
#' @export
GetHtseqSummaryStats <- function(HtseqFilePaths) {
  
  HtseqStatsPerFile <- sapply(HtseqFilePaths, function(x) { 
    table <- read.table(x, sep=":")
    counts <- table$V2
    names(counts) <- table$V1
    return(counts)
  })
  
  HtseqStatsPerFile <- as.data.frame(t(cbind(HtseqStatsPerFile)))
  rownames(HtseqStatsPerFile) <- sub("_summary_stats.txt", "", basename(HtseqFilePaths))
  return(HtseqStatsPerFile)
}
