#' Get mean read coverage from 5’ to 3’ across the most highly expressed genes. Before calling 
#' this function, select a subset of highly expressing genes. The results will be distorted 
#' if every single one of the selected genes is not highly expressed.
#'
#' Note that the genes/CDS in TopCDSDF are only considered to be suggestions of genes with high 
#' expression. For a variety of reasons (e.g. expression of certain CDS may be very rare in some 
#' genes), there may not be any expression in any of the bam files for a CDS within a gene. Those 
#' CDS will be ignored.
#'
#' @param TopCDSDF Data frame of CDS returned from select(Homo.sapiens).
#' @param PileUpPerFile List of Rsamtools pileups per file returned from GetPileUpPerGene().
#' @param NumBins Number of bins to use when parsing the counts for each gene. Default value is 100.
#' @keywords BinCounts
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A list of the median, mean, and 3rd quartile bin counts.
#' @examples
#' BamFilePaths <- dir(path="/home/ubuntu/BamFiles", pattern=".*.bam$", full.names=T)
#' 
#' TopGenesDF <- data.frame(gene = c("NM_016824", "NM_019903", "NP_001112"))
#' cdsColumns <- c("CDSCHROM", "CDSSTART", "CDSEND", "CDSSTRAND", "SYMBOL")
#' TopCDSDF <- select(Homo.sapiens, keys=as.character(TopGenesDF$gene), columns=cdsColumns, keytype="REFSEQ")
#' TopCDSDF$CDSCHROM <- sub("chr", "", TopCDSDF$CDSCHROM)
#' TopCDSDF <- TopCDSDF[which(TopCDSDF$CDSCHROM %in% c(1:22, "X", "Y")), ]
#' 
#' # BIN THE COUNTS
#' 
#' PileUpPerFile <- lapply(BamFilePaths[1:3], function(f) GetPileUpPerGene(f, TopCDSDF)) 
#' 
#' BinnedCounts <- BinCounts(TopCDSDF, PileUpPerFile, NumBins=100)
#' BinnedCountsMedian <- BinnedCounts[[1]]
#' BinnedCountsMean <- BinnedCounts[[2]]
#' BinnedCounts3rdQ <- BinnedCounts[[3]]
#' 
#' # CREATE PLOTS
#' 
#' g <- ggplot(BinnedCountsMedian, aes(x=BinID, y=RescaledCount, color=FileID)) + geom_line()
#' g <- g + ylab(label="Normalized Read Coverage") + xlab("Normalized Distance Along Transcript 5' -> 3' (%)")
#' g <- g + ggtitle(paste(format(length(unique(TopGenesDF$gene)), big.mark=","), " Genes", sep=""))
#' g
#' @export
BinCounts <- function(TopCDSDF, PileUpPerFile, NumBins=100) {
  PileUpPerGene = lapply(1:nrow(TopCDSDF), function(cdsIndex) lapply(PileUpPerFile, function(pu) pu[[cdsIndex]]));
  if (length(PileUpPerGene[[1]]) < 2) { return("ERROR: PileUpPerFile must contain at least 2 files.") }
  GenesDF = data.frame(gene=unique(TopCDSDF$REFSEQ));
  # length(which((as.character(GenesDF$gene) == as.character(TopGenesDF$gene)) == F)); # Debugging
  
  GetGeneStrand = function(ExonStrands) {
    ExonStrandTable = table(ExonStrands);
    ExonStrandTableName = names(table(ExonStrands));
    ExonStrandTableValue = table(ExonStrands);
    return(ExonStrandTableName[which.max(ExonStrandTableValue)]);
  }
  
  GenesDF$strand = unlist(lapply(GenesDF$gene, function(x) GetGeneStrand(TopCDSDF$CDSSTRAND[which(TopCDSDF$REFSEQ == x)])));
  GenesDF$strand[is.na(GenesDF$strand)] <- "*";
  
  MergeByExon = function(PileIndex) {
    # MergeByExon will create one master column of all positions "pos" found in at least one file (i.e. universe 
    # of all possible positions). Then this function will get the counts for each of these positions per file.
    NumFiles = length(PileUpPerGene[[PileIndex]]);
    p = PileUpPerGene[[PileIndex]][[1]];
    MergedDF = data.frame(pos=p$pos, count=p$count);
    for(i in 2:NumFiles) {
      p = PileUpPerGene[[PileIndex]][[i]];
      MergedDF = merge(MergedDF, data.frame(pos=p$pos, count=p$count), by="pos", all=T);
    }
    colnames(MergedDF)[seq(from=2, by=1, length.out=NumFiles)] <- c(paste("count", c(1:NumFiles), sep=""));
    MergedDF[is.na(MergedDF)] <- 0;
    MergedDF = MergedDF[order(MergedDF$pos), ];
    return(MergedDF);
  }
  
  puUnfolded = lapply(1:length(PileUpPerGene), MergeByExon);
  CountsPerCDS = unlist(lapply(puUnfolded, function(df) nrow(df)));
  # Remove CDS elements that don't have counts in any of the files.
  puCulled = puUnfolded[which(CountsPerCDS > 0)];
  TopCDSCulledDF = TopCDSDF[which(CountsPerCDS > 0), ];
  
  isValidCount = function(pos, count) {
    tryCatch(rep(pos, count),
             warning = function(w) { return(F) },
             error = function(e) { return(F) }) 
    return(T);
  }
  
  BinByExon = function(i) {
    d = puCulled[[i]];
    NumFiles = length(d);
    HistList = list();
    for(f in 2:NumFiles) {
      ValidCounts = lapply(1:nrow(d), function(j) isValidCount(d$pos[j], d[ , f][j]));
      d = d[isValidCount(d)];
      d$pos = c(1:length(d$pos)) / length(d$pos) * NumBins;
      breaks = 0:NumBins;
      flattened = unlist(lapply(1:nrow(d), function(j) rep(d$pos[j], d[ , f][j])));
      HistCounts = hist(flattened, breaks=breaks, right=F, plot=F)$counts;
      HistList[[f-1]] = HistCounts;
    }
    return(HistList);
  }
  
  puBin = lapply(1:length(puCulled), BinByExon);
  
  ReverseOrder = function(i) { 
    FilesPerGene = puBin[[i]];
    NumFiles = length(FilesPerGene);
    ReturnList = list();
    if(TopCDSCulledDF$CDSSTRAND[i]=="-") { 
      for(f in 1:NumFiles) {
        ReturnList[[f]] = rev(unlist(FilesPerGene[f])); 
      }
    } 
    else { 
      for(f in 1:NumFiles) {
        ReturnList[[f]] = unlist(FilesPerGene[f]);
      }
    }
    return(ReturnList);
  }
  
  directionalBin = lapply(1:length(puBin), ReverseOrder);
  normBin = lapply(directionalBin, function(x1) lapply(x1, function(x2) scale(x2)));
  
  BinFlatV2 = data.frame(count=unlist(normBin));
  BinFlatV2$count[is.na(BinFlatV2$count)] <- 0; # scale sometimes produces NAs
  BinFlatV2$BinID = rep(1:NumBins, length.out=nrow(BinFlatV2));
  FileIDs = 1:length(PileUpPerGene[[1]]);
  BinFlatV2$Gene = rep(TopCDSCulledDF$REFSEQ, rep(NumBins * length(FileIDs), nrow(TopCDSCulledDF)));
  BinFlatV2$FileID = rep(rep(FileIDs, rep(NumBins, length(FileIDs))), length(normBin));
  
  # For each bin of each sample, get the median/mean/3rd quartile count across all genes.
  BinPerFileMedian = ddply(BinFlatV2, .(BinID, FileID), summarise, zscore=median(count));
  BinPerFileMean = ddply(BinFlatV2, .(BinID, FileID), summarise, zscore=mean(count));
  BinPerFile3rdQ = ddply(BinFlatV2, .(BinID, FileID), summarise, zscore=as.numeric(summary(count)[5]));
  
  ScaleMinMax = function(bin) {
    ScaledMinZero = bin$zscore + abs(min(bin$zscore));
    bin$RescaledCount = ScaledMinZero / max(ScaledMinZero);
    return(bin);
  }
  
  BinPerFileMedian = rbindlist(lapply(unique(BinPerFileMedian$FileID), function(m) ScaleMinMax(BinPerFileMedian[BinPerFileMedian$FileID==m, ])));
  BinPerFileMean = rbindlist(lapply(unique(BinPerFileMean$FileID), function(m) ScaleMinMax(BinPerFileMean[BinPerFileMean$FileID==m, ])));
  BinPerFile3rdQ = rbindlist(lapply(unique(BinPerFile3rdQ$FileID), function(m) ScaleMinMax(BinPerFile3rdQ[BinPerFile3rdQ$FileID==m, ])));
  BinPerFileMedian$FileID = as.character(BinPerFileMedian$FileID);
  BinPerFileMean$FileID = as.character(BinPerFileMean$FileID);
  BinPerFile3rdQ$FileID = as.character(BinPerFile3rdQ$FileID);
  return(list(BinPerFileMedian, BinPerFileMean, BinPerFile3rdQ));
}