
source("code/utility.R")


genomeLength <- 10000
refGenome <- sample(c("A", "C", "G", "T"), size = genomeLength, replace = TRUE)
nSnp <- 50


(snpLoc <- sort(sample(1:genomeLength, size = nSnp, replace = FALSE)))


nHap <- 12
(x01 <- matrix(sample(x = 0:1, nSnp * nHap, replace = TRUE), nrow = nHap))

(snpBases <- rSnpBases(n = nSnp))

(xACTG <- convert01ToACTG(x = x01, bases = snpBases))

# THis should be a little more memory-saving, does not store the whole genome matrix
writeHaploGenomesToFasta(xACTG, snpLoc, refGenome, ploidy = 2, outDir = "haploGenomes")
writeHaploGenomesToFasta(xACTG, snpLoc, refGenome, ploidy = 4, outDir = "haploGenomes4")
