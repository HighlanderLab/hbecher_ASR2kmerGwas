## GG's functions ##

# Converting AlphaSimR SNP haplotype matrix (0/1) to A/C/T/G and full genome.
# The aim is to be able to generate k-mers from this genome data using a tool like
# https://sourceforge.net/projects/kanalyze.

# TODO: This is a prototype implementation and is likely inefficient - identify
#       bottlenecks and optimise.

genomeLength <- 10
refGenome <- sample(c("A", "C", "G", "T"), size = genomeLength, replace = TRUE)
nSnp <- 3
pos <- c(2, 5, 7)
snpLoc <- sort(sample(1:genomeLength, size = nSnp, replace = FALSE))
nHap <- 4
x01 <- matrix(sample(x = 0:1, nSnp * nHap, replace = TRUE), nrow = nHap)

rSnpBases <- function(n) {
  #' Generate a random set of n bi-allelic SNPs bases with transition mutations
  #' @n number of SNPs
  #' @return 2-column matrix holding ancestral alleles in column 1 and derived alleles in column 2
  #' @examples
  #' rSnpBases(n = 20)
  # Only bi-allelic SNP, transition mutations, and showing one strand only
  baseMutations <- matrix( # fmt: skip
    data = c("A", "G",
             "C", "T",
             "G", "A",
             "T", "C"),
    nrow = 4, ncol = 2, byrow = TRUE) # fmt: skip
  ancestralBaseInt <- sample(x = 1:4, size = n, replace = TRUE)
  ancestralBaseChar <- baseMutations[ancestralBaseInt, 1]
  derivedBaseChar <- baseMutations[ancestralBaseInt, 2]
  bases <- matrix(
    data = c(ancestralBaseChar, derivedBaseChar),
    ncol = 2,
    byrow = FALSE
  )
  return(bases)
}


convert01ToACTG <- function(x, bases) {
  #' Convert a matrix of haplotypes encoded as 0/1 to A/C/T/G assuming bi-allelic SNPs
  #' @x numeric matrix of 0 and 1
  #' @bases 2-column matrix holding ancestral alleles in column 1 and derived alleles in column 2
  #' @return character matrix of A, C, T, and G
  #' @examples
  #' x <- matrix(sample(x = 0:1, 20, replace = TRUE), nrow = 4)
  #' convert01ToACTG(x, bases = rSnpBases(n = 20))
  nLoc <- ncol(x)
  ret <- matrix(data = "", nrow = nrow(x), ncol = nLoc)
  for (loc in 1:nLoc) {
    ret[, loc] = bases[loc, ][x[, loc] + 1]
  }
  return(ret)
}


expandSnpToGenome <- function(x, loc, genome) {
  #' Expand a SNP haplotype matrix to a full genome sequence matrix
  #' @x character matrix of A/C/T/G bases
  #' @loc integer vector indicating SNP locations in the genome
  #' @genome character vector of the full genome sequence
  #' @return character matrix of the full genome sequence
  #' @examples
  #' x <- matrix(sample(c("A", "C", "G", "T"), size = 20, replace = TRUE), nrow = 4)
  #' loc <- sort(sample(1:100, size = 20, replace = FALSE))
  #' genome <- sample(c("A", "C", "G", "T"), size = 100, replace = TRUE)
  #' expandSegSitesToGenome(x, loc, genome)
  nInd <- nrow(x)
  genomeLength <- length(genome)
  ret <- matrix(data = "", nrow = nInd, ncol = genomeLength)
  for (ind in 1:nInd) {
    ret[ind, ] <- genome
    ret[ind, loc] <- x[ind, ]
  }
  return(ret)
}




# function to generate a genome reference, vectore of character
makeRefGenome <- function(length) {
  #' Generate a random reference genome sequence
  #' @length integer length of the genome
  #' @return character vector of the full genome sequence
  #' @examples
  #' makeRefGenome(length = 100)
  genome <- sample(c("A", "C", "G", "T"), size = length, replace = TRUE)
  return(genome)
}

# function to write out haplogenomes based on reference and snps, to fasta files,
# one individual per file, this should be ploidy-aware
writeHaploGenomesToFasta <- function(
    x,
    loc,
    genome,
    prefix = "individual",
    outDir = "haploGenomes",
    ploidy = 2
) {
  #' Write out haplogenomes to fasta files, one individual per file
  #' @x character matrix of A/C/T/G bases
  #' @loc integer vector indicating SNP locations in the genome
  #' @genome character vector of the full genome sequence
  #' @prefix character prefix for the output fasta files
  #' @outDir character output directory for the fasta files
  #' @ploidy integer ploidy of the individuals
  #' @return NULL
  #' @examples
  #' x <- matrix(sample(c("A", "C", "G", "T"), size = 20, replace = TRUE), nrow = 4)
  #' loc <- sort(sample(1:100, size = 20, replace = FALSE))
  #' genome <- sample(c("A", "C", "G", "T"), size = 100, replace = TRUE)
  #' writeHaploGenomesToFasta(x, loc, genome, ploidy = 2)
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  nInd <- nrow(x) / ploidy
  for (ind in 1:nInd) {
    haploGenomes <- expandSnpToGenome(
      x = x[((ind - 1) * ploidy + 1):(ind * ploidy), ],
      loc = loc,
      genome = genome
    )
    fastaLines <- c()
    for (haplo in 1:ploidy) {
      fastaLines <- c(fastaLines, paste0(">", prefix, ind, "_haplo", haplo))
      fastaLines <- c(fastaLines, paste(haploGenomes[haplo, ], collapse = ""))
    }
    writeLines(fastaLines, con = file.path(outDir, paste0(prefix, ind, ".fasta")))
  }
  return(NULL)
}
