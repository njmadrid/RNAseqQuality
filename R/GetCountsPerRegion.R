#' Count the number of reads that are in a given genomic region.
#' 
#' @param BamFilePath location of aligned bam file. Corresponding bai index file must be in the same directory.
#' @param param ScanBamParam data structure.
#' @keywords GetCountsPerRegion
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return An integer value for the number of reads that are in a given genomic region.
#' @examples
#' # Get file paths of aligned bam files. The corresponding bai index file must also be in the same directory.
#' BamFilePaths <- dir(path="/home/BamFiles", pattern=".*.bam$", full.names=T)
#' rRNAparamRegions <- GRanges(seqnames="GL000220.1", IRanges(1, 161802))
#' 
#' # Note that any region(s) can be used in GetCountsPerRegion, not just rRNA regions. 
#' rRNAparam <- ScanBamParam(which=rRNAparamRegions)
#' 
#' # Secondary alignments will only be counted once
#' rRNATotalCounts <- GetCountsPerRegion(BamFilePaths[1], rRNAparam)
#' @export
GetCountsPerRegion <- function(BamFilePath, param) {
  df <- readGAlignments(BamFilePath, param=param, use.names=T)
  TotalCount = length(unique(names(df))) # don't count secondary alignments
  return(TotalCount)
}