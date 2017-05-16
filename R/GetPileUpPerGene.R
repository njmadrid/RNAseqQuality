#' Get read depths per gene as defined by the RSamtools pileup function.
#'
#' @param BamFilePath Location of aligned bam file. Corresponding bai index file must be in the same directory.
#' @param TopCDSDF Data frame of CDS returned from select(Homo.sapiens).
#' @param max_depth Parameter of the same name passed to Rsamtools pileup. Default value is 10000.
#' @keywords GetPileUpPerGene
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame created via Rsamtools:pileup for each CDS region in TopCDSDF.
#' @examples
#' BamFilePaths <- dir(path="/home/ubuntu/BamFiles", pattern=".*.bam$", full.names=T)
#' TopGenesDF <- data.frame(gene = c("NM_016824", "NM_019903", "NP_001112"))
#' cdsColumns <- c("CDSCHROM", "CDSSTART", "CDSEND", "CDSSTRAND", "SYMBOL")
#' TopCDSDF <- select(Homo.sapiens, keys=as.character(TopGenesDF$gene), columns=cdsColumns, keytype="REFSEQ")
#' TopCDSDF$CDSCHROM <- sub("chr", "", TopCDSDF$CDSCHROM)
#' TopCDSDF <- TopCDSDF[which(TopCDSDF$CDSCHROM %in% c(1:22, "X", "Y")), ]
#' PileUpPerFile <- lapply(BamFilePaths[1:3], function(f) GetPileUpPerGene(f, TopCDSDF)) 
#' @export
GetPileUpPerGene = function(BamFilePath, TopCDSDF, max_depth=10000) {
  TopCDSGR = GRanges(seqnames=as.character(TopCDSDF$CDSCHROM), ranges=IRanges(start=TopCDSDF$CDSSTART, end=TopCDSDF$CDSEND));
  
  GetPileUpPerExon = function(which) {
    pu = pileup(BamFilePath,
                scanBamParam=ScanBamParam(which=which),
                pileupParam=PileupParam(max_depth=max_depth,
                                        distinguish_strands=F,
                                        distinguish_nucleotides=F,
                                        include_deletions=F,
                                        include_insertions=F));
    return(pu);
  }
  
  puList = lapply(TopCDSGR, GetPileUpPerExon);
  
  return(puList);
}