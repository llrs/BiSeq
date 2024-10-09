.annotateGRanges <- function(object, regions, name, regionInfo){

  regions <- unique(regions)
  MetaCol <- ncol(mcols(object)) + 1
  if(missing(regionInfo)){
    ind <- overlapsAny(object, regions)
    mcols(object)[, MetaCol] <- FALSE
    mcols(object)[, MetaCol][ind] <- TRUE
    colnames(mcols(object))[MetaCol] <- name
    return(object)
  }
  
  if(is(mcols(regions)[,regionInfo], "factor")){
    mcols(regions)[, regionInfo] <- as.character(mcols(regions)[, regionInfo])
  }
  overl <- findOverlaps(query=object, subject=regions)
  matches <- as.data.frame(overl)
  
  helper <- rep(NA, length(object))
  if(nrow(matches) > 0){
    tab <- as.data.frame(table(matches$queryHits)) # how many regions map to region in object
    
    ids <- mcols(regions)[, regionInfo]
    if(is(ids, "CompressedCharacterList")){
      ids.l <- sapply(ids, length)
      ids.char <- character(length=length(ids))
      ids.l.1 <- which(ids.l == 1)
      ids.l.n <- which(ids.l > 1)
      if(length(ids.l.1) > 1){
        ids.char[ids.l.1] <- unlist(ids[ids.l.1])
      }
      if(length(ids.l.n) > 1){
        ids.char[ids.l.n] <- sapply(ids[ids.l.n], function(x) paste(x, collapse = ","))
      }
      ids <- ids.char
    } 
    
    ind.many <- as.numeric(as.character(tab$Var1[tab$Freq > 1]))
    matches.one <- matches[!is.element(matches$queryHits, ind.many), ]
    
    helper[matches.one$queryHits] <- ids[matches.one$subjectHits]
    
    # If there aren't multiple mapped positions
    if (length(ind.many) == 0L) {
      mcols(object)[, MetaCol] <- helper
      colnames(mcols(object))[MetaCol] <- name
      return(object)
    }
    
    for(i in ind.many){
      ind.reg <- matches$subjectHits[matches$queryHits == i]
      names.reg <- ids[ind.reg]
      names.reg <- sort(unique(names.reg))
      ids.i <- paste(names.reg, collapse=",")
      helper[i] <- ids.i
    }
  }
  mcols(object)[, MetaCol] <- helper
  colnames(mcols(object))[MetaCol] <- name
  object
}


setMethod("annotateGRanges",
    signature=c(object="GRanges", regions="GRanges", name="character", regionInfo="character"),
    .annotateGRanges)

setMethod("annotateGRanges",
    signature=c(object="GRanges", regions="GRanges", name="character", regionInfo="integer"),
    .annotateGRanges)

setMethod("annotateGRanges",
    signature=c(object="GRanges", regions="GRanges", name="character", regionInfo="missing"),
    .annotateGRanges)
