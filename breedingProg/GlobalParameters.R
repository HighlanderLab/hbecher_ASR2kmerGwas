# Global Parameters

# ---- Number of simulation replications and breeding cycles ----
nReps = 1 # Number of simulation replicates
nBurnin = 10 # Number of years in burnin phase
nFuture = 3 # Number of years in future phase
nCycles = nBurnin + nFuture
startTP = 10 # Year to start collecting genomic data

# ---- Genome simulation ----
quickGenomeSim = FALSE # it doesn't make much of a difference
#nChr = 16
nChr = 3
#nQtl = 1000 # Number of QTL per chromosome
nQtl = 10 # Number of QTL per chromosome
#nSnp = c(256, 170000) # Number of SNP per chromosome
nSnp = c(100) # Number of SNP per chromosome
ploidy=2

# ---- Initial parents mean and variance ----
initMeanG = 1 # Initial mean of genetic values
initVarG = 1 # Initial additive genetic variance
initRelAA = 0.2 # relative additive-by-additive variance compared to initVarG
#                 (in a diploid organism with allele frequency 0.5)
initVarEnv = 1e-6 # Virtually zero for consistency with 2-Part paper
initVarGE = 0 # Genotype-by-Environment (GxE) variance (setting to 0 to keep things simple)
varE = c(4, 4) # Yield trial error variance, bushels per acre
# Relates to error variance for an entry mean

# ---- Breeding program details ----
nParents = 50 # Number of parents to start a breeding cycle
nCrosses = 100 # Number of crosses per year
nDH = 100 # DH lines produced per cross
famMax = 10 # The maximum number of DH lines per cross to enter PYT
nPYT = 1000 # Entries per preliminary yield trial (TRAINING DATA)
nAYT = 100 # Entries per advanced yield trial (TRAINING DATA)
nEYT = 10 # Entries per elite yield trial (TRAINING DATA)

# Effective replication of yield trials
repHDRW = 4 / 9 # h2 = 0.1
repPYT = 1 # h2 = 0.2
repAYT = 4 # h2 = 0.5
repEYT = 8 # h2 = 0.7
