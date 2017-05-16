#' Get overlaps using GenomicRanges findOverlaps for a bedgraph coverage file (*.bedgraph.gz) compared to RefSeq defined gene regions.
#'
#' @param bedgraph GRanges of bedgraph coverage file (*.bedgraph.gz) created via bedtools genomecov -bga -split.
#' @param refseq GRanges of refseq defined gene regions.
#' @keywords GetBedGraphRefSeqOverlaps
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return Results are returned as a Hits object.
#' @examples
#' refseq <- import.gff("RefSeq_hg19_exons_nodups_021114g1k.gff")
#' refseq$gene <- sub("_exon_[0-9]+", "", refseq$gene_id)
#' refseq$width <- width(refseq)
#' txlength <- tapply(refseq$width, refseq$gene, sum)

#' # Get file paths of aligned bedgraph files
#' # Coverage calculated w/ bedtools genomecov -bga -split
#' BedGraphFiles <- dir(path="/home/bedgraphs", full.names=T)
#' BedgraphList <- lapply(BedGraphFiles[1:3], function(bg) import.bedGraph(bg))
#' BedGraphRefSeqOverlapsList <- lapply(BedgraphList, function(bg) GetBedGraphRefSeqOverlaps(bg, refseq))
#' @export
GetBedGraphRefSeqOverlaps <- function(bedgraph, refseq) {
  # this probably isn't necessary, but some GRanges operations whine otherwise
  seqlevels(refseq) <- seqlevels(bedgraph)
  # get the coverage info overlapping exons
  BedGraphRefSeqOverlaps <- subsetByOverlaps(bedgraph, refseq)
  # track the gene names back on those ranges
  BedGraphRefSeqOverlaps$gene <- mcols(refseq)[findOverlaps(BedGraphRefSeqOverlaps, refseq, select="first"),]$gene
  # we have the over lapping ranges, but they haven't been trimmed to fit the exons
  refhit <- findOverlaps(BedGraphRefSeqOverlaps, refseq, select="first")
  start(BedGraphRefSeqOverlaps) <- ifelse(start(BedGraphRefSeqOverlaps) >= start(refseq[refhit]), start(BedGraphRefSeqOverlaps), start(refseq[refhit]))
  end(BedGraphRefSeqOverlaps) <- ifelse(end(BedGraphRefSeqOverlaps) <= end(refseq[refhit]), end(BedGraphRefSeqOverlaps), end(refseq[refhit]))
  return(BedGraphRefSeqOverlaps)
}
