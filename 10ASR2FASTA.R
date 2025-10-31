

rm(list = ls())
library(here)
library(AlphaSimR)


breedingProgDir <- "breedingProg"


# Parameters --------------------------------------------------------------


# Reference genome/ancestral alleles --------------------------------------

chrLens <- c(10000,10000, 20000) # in nt

refList = lapply(chrLens, \(x) sample(c("A","C","G","T"), size=x, replace=TRUE))


# Run breeding simulation ------------------------------------------------

source(paste(here(), breedingProgDir, "00RUNME.R", sep="/"))



# Generate variant data (haplotype 'matrix') ---------------------------------------------------



# genet map is a list of vectors
gMap <- SP$genMap
# number of 'variant' sites per chromosome
varSitesPerChr <- sapply(gMap, length)

cMap <- do.call(c, gMap) # concatenated map (relies on element names)

trait1lpc <- SP$traits[[1]]@lociPerChr

trait1ll <- SP$traits[[1]]@lociLoc

sum(trait1lpc) == length(trait1ll)

# cumulative sums of number of QTLs per chromosome
cs <- c(0, cumsum(trait1lpc))

# cumulative sums of number of var sites per chromosome
csChrom <- c(0, cumsum(varSitesPerChr))

# a list of vectors, one for each chromosome, with the indices of the QTLs (among all variant sites)
trait1chromPosList <- lapply(1:(length(cs)-1), function(x) {
  trait1ll[(cs[x]+1):cs[x+1]]+csChrom[x]
})

SP$traits[[1]]@addEff
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


t1data <- getTraitQtlData(1, SimParam=SP)

# gpl - list of genetic positions (one vector per chromosome)
# chrl - vector of chromosome lenghts
# mfl - (optional) a list of mapping functions (gneet pos to nt pos)
genetPos2ntPos <- function(gpl, chrl, mfl){
  if(!missing(mfl)) stop("Mapping functions are not implemented yet.")
  stopifnot(length(gpl) == length(chrl))
  intPosList <- mapply(function(pos, len) {
    round(pos/ max(pos) * len)
  }, gpl, chrl, SIMPLIFY = FALSE)
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

ntPos <- genetPos2ntPos(SP$genMap, chrLens)

# none duplicated?
sapply(ntPos, \(x) any(duplicated(x)))



# TODO check what population(s) to use here
haplo <- pullSegSiteHaplo(AYT, simParam = SP)
dim(haplo)

haplo2list <- function(haplo, varSitesPerChr){
  lst <- list()
  cs <- c(0, cumsum(varSitesPerChr))
  for(i in 1:length(varSitesPerChr)){
    lst[[i]] <- haplo[ , (cs[i]+1):cs[i+1] ]
  }
  return(lst)
}

haploList <- haplo2list(haplo, varSitesPerChr)
lapply(haploList, dim)




# decide which variant sites are indels
prIns <- 0.1 # probability of an insertion at a variant site
makeInsertionBoolList <- function(vcpc, prIns=prIns){
  lapply(1:length(vcpc), \(x){
  sample(c(T, F), prob = c(prIns, 1-prIns), size=vcpc[x], replace = TRUE)
})
}

insBoolList <- makeInsertionBoolList(varSitesPerChr, prIns=0.1)


# a function to create insertion alleles of a given length
makeInsertionAllele <- function(len){
  paste0(sample(c("A","C","G","T"), size=len, replace = TRUE), collapse = "")
}
#makeInsertionAllele(30)


# a function to make alternate alleles for each variant site, taking into
#  account the reference allele and whether it's going to be an insertion
makeAltAlleles <- function(ref, pos, isInsertion){
  altAlleles <- list()  
  for(i in 1:length(pos)){
    if(isInsertion[i]) {
      altAlleles[i] <- makeInsertionAllele(30)
    } else {
      altAlleles[i] <- sample(setdiff(c("A", "C", "G", "T"), ref[pos[i]]  ), 1)
    }
    
  }
  return(altAlleles)
}

# TODO GO ON HERE

# TODO, make work for multiple chromosomes
# make some alt alleles to be used below
aa <- makeAltAlleles(ref=ref, pos=pos, isInsertion = isInsertion)

# one haplotype, a string of length >= LL (due to insertions) 
makeFastaString <- function(ref, hap, pos, aa, hapIndex){
  vec <- ref # vector as long as the ref genome
  
  # indexing over number of variant sites
  #  replace allele if corresponding index in the respective haplotype is 1
  #  (i.e., if derived state)
  for(i in 1:length(pos)){
    if(hap[hapIndex,i] == 1) vec[pos[i]] <- aa[i]
  }
  return(paste0(paste0(vec, collapse=""), "\n"))
}


# function to write out a FASTA file of the deired ploidy
#  ploidy is set by the length of idx (a vector)
makeFasta <- function(fname, idx=c(1,2), ref, hap, aa, pos){
  f <- file(fname, "w")
  cat("", file=fname)
  for(i in 1:length(idx)){
    print(i)
  
    cat(paste0(">haplotype", idx[i], "\n"), file=fname, append = TRUE)
    cat(makeFastaString(ref=ref, hap=hap, pos=pos, aa=aa, hapIndex=i), file=fname, append=TRUE )
  }
  close(f)
}



# Write out files ---------------------------------------------------------

#create FASTA files to be used as the input for kmc
makeFasta("my12.fa", idx=c(1,2), ref=ref, hap=haplo, aa, pos)
makeFasta("my34.fa", idx=c(1,2), ref=ref, hap=haplo, aa, pos)

