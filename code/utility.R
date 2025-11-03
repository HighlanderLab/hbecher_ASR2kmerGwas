## GG's functions ##

# Converting AlphaSimR SNP haplotype matrix (0/1) to A/C/T/G and full genome.
# The aim is to be able to generate k-mers from this genome data using a tool like
# https://sourceforge.net/projects/kanalyze.

# TODO: This is a prototype implementation and is likely inefficient - identify
#       bottlenecks and optimise.


# QTL indices and additive effects per chromosome (two lists))
getTraitQtlData <- function(trt, SimParam=SP){
  traitlpc <- SimParam$traits[[trt]]@lociPerChr
  traitll <- SimParam$traits[[trt]]@lociLoc
  
  cs <- c(0, cumsum(traitlpc))
  
  traitchromPosList <- lapply(1:(length(cs)-1), function(x) {
    
    traitll[(cs[x]+1):cs[x+1]]})
  
  posList <- mapply(function(x,y) {
    x[y]
  }, SP$genMap, traitchromPosList, SIMPLIFY = FALSE)
  
  addEffList <- lapply(1:(length(cs)-1), function(x) {
    SP$traits[[trt]]@addEff[(cs[x]+1):cs[x+1]]
  })
  
  return(list(pos=posList, addEff=addEffList))
}



# gpl - list of genetic positions (one vector per chromosome)
# chrl - vector of chromosome lenghts
# mfl - (optional) a list of mapping functions (gneet pos to nt pos)
genetPos2ntPos <- function(gpl, chrl, mfl){
  if(!missing(mfl)) stop("Mapping functions are not implemented yet.")
  stopifnot(length(gpl) == length(chrl))
  intPosList <- mapply(function(pos, len) {
    round(pos/ max(pos) * len)+1
  }, gpl, chrl-1, SIMPLIFY = FALSE)
  # space out positions to avoid duplicates
  intPosList <- lapply(intPosList, function(x) {
    dupIdx <- which(duplicated(x))
    while(length(dupIdx) > 0){
      x[dupIdx] <- x[dupIdx] + 1
      dupIdx <- which(duplicated(x))
    }
    return(x)
  })
  return(intPosList)
}



# Convert haplotype matrix to list of matrices (one per chromosome)
haplo2list <- function(haplo, varSitesPerChr){
  lst <- list()
  cs <- c(0, cumsum(varSitesPerChr))
  for(i in 1:length(varSitesPerChr)){
    lst[[i]] <- haplo[ , (cs[i]+1):cs[i+1] ]
  }
  return(lst)
}



# logical indicator vectors for which variant sites are insertions
makeInsertionBoolList <- function(vcpc, prIns=prIns){
  lapply(1:length(vcpc), \(x){
    sample(c(T, F), prob = c(prIns, 1-prIns), size=vcpc[x], replace = TRUE)
  })
}



# a function to create insertion alleles of a given length
makeInsertionAllele <- function(len){
  paste0(sample(c("A","C","G","T"), size=len, replace = TRUE), collapse = "")
}
#makeInsertionAllele(30)


# a function to make alternate alleles for each variant site, taking into
#  account the reference allele and whether it's going to be an insertion
# ref, po, and isInsertion are lists of vectors
makeAltAlleles <- function(ref, pos, isInsertion){
  altAlleles <- list()  
  
  for(i in 1:length(pos)){
    currList <- character()
    for(j in 1:length(pos[[i]])){

      if(isInsertion[[i]][j]) {
        currList[j] <- makeInsertionAllele(30)
      } else {
        currList[j] <- sample(setdiff(c("A", "C", "G", "T"), ref[[i]][pos[[i]][j]]  ), 1)
      }
#      print(currList[j])
    }
    altAlleles[[i]] <- currList
  }
  return(altAlleles)
}


# a function to consolidate locus and QTL infomration
makeVarLocList <- function(tDat, posList, isIns){
  locDf <- data.frame(id=do.call(c, lapply(posList, \(x) names(x))),
                      chr=rep(1:length(posList), times=sapply(posList, length)),
                      loc=unname(do.call(c, posList)),
                      eff=0,
                      isIns=unname(do.call(c, isIns))
  )
  for(i in 1:length(tDat$pos)){
    for(j in 1:length(tDat$pos[[i]])){
      locDf$eff[locDf$id == names(tDat$pos[[i]][j])] <-
        tDat$addEff[[i]][j]
    }
  }
  return(locDf)
}




# one haplotype for a single chromosome, a string
makeFastaString <- function(ref, hapArr, pos, aa, hapIndex){
  vec <- ref # vector as long as the ref genome
  
  # indexing over number of variant sites
  #  replace allele if corresponding index in the respective haplotype is 1
  #  (i.e., if derived state)
  for(i in 1:length(pos)){
    if(hapArr[hapIndex,i] == 1) vec[pos[i]] <- aa[i]
  }
  return(paste0(paste0(vec, collapse=""), "\n"))
}


# one haplotype for each chromosome, a list of strings
makeFastaStringList <- function(refL, hapArrL, posL, aaL, hapIndex){
  fastaStrL <- list()
  for(chr in 1:length(refL)){
    fastaStrL[[chr]] <- makeFastaString(ref=refL[[chr]],
                                        hapArr=hapArrL[[chr]],
                                        pos=posL[[chr]],
                                        aa=aaL[[chr]],
                                        hapIndex=hapIndex)
  }
  return(fastaStrL)
}

# function to write out a FASTA file of the desired ploidy level
#  ploidy is set by the length of the vector idx
#  uses makeFastaStringList to generate afst strings
makeFasta <- function(fname, idx=c(1,2), ref, hap, aa, pos, varSitesPerChrom, indName){
  f <- file(fname, "w")
  for(hpt in 1:length(idx)){
    fsl <- makeFastaStringList(refL=ref,
                                 hapArrL=haplo2list(hap, varSitesPerChrom),
                                 posL=pos,
                                 aaL=aa,
                                 hapIndex=idx[hpt])
    for(chr in 1:length(fsl)){
      cat(paste0(">Ind_", indName, "_hap_", hpt, "_chr_", chr, "\n"),
          file=fname, append=TRUE)
      cat(fsl[[chr]], file=fname, append=TRUE)
    }
   }
  close(f)
}


