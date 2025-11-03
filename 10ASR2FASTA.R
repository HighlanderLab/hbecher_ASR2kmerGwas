

rm(list = ls())
library(here)
library(AlphaSimR)

source(here("code/utility.R"))

breedingProgDir <- "breedingProg"
fastaDir <- "fasta"

if (!dir.exists(here(fastaDir))) {
  dir.create(here(fastaDir), recursive = TRUE)
}



# Parameters --------------------------------------------------------------
prIns <- 0.1 # probability of an insertion at a variant site

# chromosome lengths
chrLens <- c(30000,30000,30000) # in nt

# reference alleles
refList = lapply(chrLens, \(x) sample(c("A","C","G","T"), size=x, replace=TRUE))


# Run breeding simulation ------------------------------------------------

source(here(breedingProgDir, "00RUNME.R"))


# the simulation writes out 4 haplotype files per future generation (TXT.GZ) and
# 4 pheno files (TXT)


# we need code to turn these haplotype files into FASTA files
# this will depend on the reference, alternative alleles, marker locations, and
# haplotype data



# collect site information ---------------------------------------------------

# genet map is a list of vectors
gMap <- SP$genMap
# number of 'variant' sites per chromosome
varSitesPerChr <- sapply(gMap, length)

cMap <- do.call(c, gMap) # concatenated map (relies on element names)

trait1lpc <- SP$traits[[1]]@lociPerChr

trait1ll <- SP$traits[[1]]@lociLoc

sum(trait1lpc) == length(trait1ll)






SP$traits[[1]]@addEff


t1data <- getTraitQtlData(1, SimParam=SP)

ntPos <- genetPos2ntPos(SP$genMap, chrLens)


# none duplicated?
sapply(ntPos, \(x) any(duplicated(x)))


# decide which variant sites are indels
insBoolList <- makeInsertionBoolList(varSitesPerChr, prIns=0.1)



# make some alt alleles to be used below
aa <- makeAltAlleles(ref=refList, pos=ntPos, isInsertion = insBoolList)


# Table of var sites and effect sizes -------------------------------------

tDat <- getTraitQtlData(1, SimParam=SP)
ntList <-  genetPos2ntPos(SP$genMap, chrLens)

locusTable <- makeVarLocList(tDat, ntList, insBoolList)
write.table(locusTable, file="AllVarLoci.txt", row.names = FALSE, quote=FALSE)
# Consider adding k-mer sequences (ref and alt) to the table,
#  could be comma-separated strings
#  this will depend on the k-mer length




# Consolidate pheno files produced ----------------------------------------

phenoFiles <- dir(pattern="phenoPredict")
phenoData <- do.call(rbind, lapply(phenoFiles, function(f){
  read.table(f, header=TRUE)
}))
write.table(phenoData, file="AllPhenoData.txt", sep="\t", row.names=FALSE, quote=FALSE)
head(phenoData)

# same thin for gv files
gvFiles <- dir(pattern="gvPredict")
gvData <- do.call(rbind, lapply(gvFiles, function(f){
  read.table(f, header=TRUE)
}))
write.table(gvData, file="AllGvData.txt", sep="\t", row.names=FALSE, quote=FALSE)


#pairs(cbind(gvData, phenoData)[,-c(1, 3, 5, 7, 8, 10)])

# Create FASTAs -----------------------------------------------------------

haploFiles <- dir(pattern="haplotypes")

# if these are too many files, focus on a subset only
haploFiles <- dir(pattern="haplotypes*PYT")

# loop over haplo files
for(hf in haploFiles){
  print(paste0("Working on file ", hf, "..."))
  
  hapMat <- as.matrix(read.table(gzfile(hf), header=FALSE, row.names=1))
  
  
  for(ind in 1:(nrow(hapMat)/ploidy)){
    indNam = strsplit(row.names(hapMat)[(ind-1)*ploidy +1], split="_")[[1]][1]
    
    fname <- here(fastaDir, paste0("Ind", indNam, ".fa"))
    print(paste0(" Writing haplotype FASTA to ", fname, "..."))
    makeFasta(fname,
              idx=c((ind-1)*ploidy +1,(ind-1)*ploidy +2),
              ref=refList,
              hap=hapMat,
              aa,
              ntList,
              varSitesPerChr=varSitesPerChr,
              indName=indNam)
    
  }
  
  
}
