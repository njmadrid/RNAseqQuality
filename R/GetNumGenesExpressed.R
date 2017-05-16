#' Get the total number of genes expressed per gtf file (generated via Cufflinks).
#'
#' @param gtfDFs data table structure containing Cufflinks results.
#' @param covThreshold only genes will be count that have "cov" values above this threshold. Default value is 1.0.
#' @keywords GetNumGenesExpressed
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return An integer for the total number of genes expressed per gtf file (generated via Cufflinks).
#' @examples
#' gtffiles <- c("transcripts_sorted.gtf")
#' FileCons <- mclapply(gtffiles, function(FileName) file(FileName), mc.cores=8)
#' gtfs <- mclapply(FileCons, function(fc) import.gff(fc, format="gtf"), mc.cores=8)
#' for(i in length(FileCons)) { close(FileCons[[i]]) }
#' gtfDFs <- mclapply(gtfs, function(gtf) as.data.frame(gtf[gtf$type=="transcript"]), mc.cores=8)
#' names(gtfDFs) <- basename(gtffiles)
#' NumGenesPerFile <- GetNumGenesExpressed(gtfDFs, 1.0)
#' @export
GetNumGenesExpressed <- function(gtfDFs, covThreshold=1.0) {
  fixedgtfs <- foreach (i=iter(gtfDFs), n=iter(names(gtfDFs))) %do% cbind(FileID = n, i) # Add the FileID column
  cufftbl <- rbindlist(fixedgtfs)
  isoforms.expressed <- tapply(cufftbl$cov, cufftbl$FileID, function(x) sum(x>=covThreshold))
  isAcceptable <- tapply(cufftbl$cov, paste(cufftbl$FileID, cufftbl$gene_id, sep="_"), function(x) any(x>=covThreshold))
  genes.expressed <- tapply(isAcceptable, sub("_ENSG[0-9]+$","", names(isAcceptable)), sum)
  return(genes.expressed)
}
