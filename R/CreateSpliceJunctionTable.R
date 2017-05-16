#' Create a table of unique splice junction expression per tab file (*.out.tab) generated using during
#' some alignment processes (e.g. using the STAR aligner).
#'
#' @param tabfiles directory path of *.out.tab files generated during alignment.
#' @keywords CreateSpliceJunctionTable
#' @author 
#' Nathaniel J. Madrid, Jason Byars
#' @return A data.frame of the unique splice junction expression per tab file (*.out.tab).
#' @examples
#' tabfiles <- dir(path="tabs/", pattern=".*.out.tab$", full.names=T)
#' JunctionTables <- CreateSpliceJunctionTable(tabfiles)
#' @export
CreateSpliceJunctionTable <- function(tabfiles) {
  # create a merged splice junction table
  TabTable <- lapply(tabfiles, read.table, header=F, stringsAsFactors=F)
  names(TabTable) <- sub("^rna_","",sub("_SJ.out.tab$", "", basename(tabfiles)))
  junctiontbl <- bind_rows(TabTable, .id = "sample")
  colnames(junctiontbl) <- c("sample", "chromosome","first.intron.base","last.intron.base", "strand", "intron.motif", 
                             "annotated", "unique.mapping","multi.mapping","max.overhang")
  
  unique.mapping <- junctiontbl %>% select(sample,chromosome, first.intron.base, last.intron.base, 
                                           strand, intron.motif, annotated, unique.mapping) %>% 
    spread(key=sample, value=unique.mapping, fill=0)
  
  multi.mapping <- junctiontbl %>% select(sample,chromosome, first.intron.base, last.intron.base, 
                                          strand, intron.motif, annotated, multi.mapping) %>% 
    spread(key=sample, value=multi.mapping, fill=0)
  
  all.mapping <- junctiontbl %>% mutate(all.mapping = unique.mapping + multi.mapping) %>%
    select(sample,chromosome, first.intron.base, last.intron.base, 
           strand, intron.motif, annotated, all.mapping) %>% 
    spread(key=sample, value=all.mapping, fill=0)
  
  # how many novel vs annotated junctions
  table(all.mapping$annotated) # 670302 vs 265411
  # granted a good portion of the novel will be garbage
  
  
  # get a hint what a lower threshold should be. 5 reads to count an event appears reasonable
  table(junctiontbl$unique.mapping)[1:20]
  table(junctiontbl$multi.mapping)[1:20]
  table(junctiontbl$unique.mapping+junctiontbl$multi.mapping)[1:20]
  
  # junctions observed
  total.junctions <- apply(all.mapping[,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))
  annotated.junctions <- apply(all.mapping[all.mapping$annotated==1,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))
  novel.junctions <- apply(all.mapping[all.mapping$annotated==0,7:dim(all.mapping)[2]], 2, function(x) sum(x>=5))
  
  df <- data.frame(total.junctions, annotated.junctions, novel.junctions)
  
  return(df)
}