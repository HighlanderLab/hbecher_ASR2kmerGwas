
setwd("~/git_repos/hbecher_ASR2kmerGwas/")

library(AlphaSimR)



# Parameters --------------------------------------------------------------


mu=1e-7
LL=1e5
prIns <- 0.1 # proportion of mutation that are insertions


# Reference genome/ancestral alleles --------------------------------------

ref = sample(c("A","C","G","T"), size=LL, replace=TRUE)


# Generate variant data (haplotype 'matrix') ---------------------------------------------------
founders <- runMacs2(nInd=1,
                     nChr=1,
                     segSites=NULL,
                     Ne=200000,
                     bp=LL,
                     genLen = LL/1e8,
                     ploidy=2,
                     mutRate=mu)

mySP <- SimParam$new(founderPop = founders)
pop0 <- newPop(rawPop = founders,
               simParam = mySP)
haplo <- pullSegSiteHaplo(pop0, simParam = mySP)


# Prepare for conversion --------------------------------------------------
# Convert map locations to integer nucleotide positions, making sure there are no duplicated sites created during conversion!
stopifnot(all(!duplicated(as.integer(mySP$genMap$"1"*1e8)))); pos <- as.integer(mySP$genMap$"1"*1e8)

# decide which variant sites are indels
isInsertion <- sample(c(T, F), prob = c(prIns, 1-prIns), size=length(pos), replace = TRUE)

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

