# Create founders

# # Generate initial haplotypes
# if (quickGenomeSim) {
#   founderPop = quickHaplo(
#     nInd = nParents,
#     nChr = nChr,
#     segSites = nQtl + max(nSnp),
#     inbred = TRUE
#   )
# } else {
#   founderPop = runMacs(
#     nInd = nParents,
#     nChr = nChr,
#     segSites = nQtl + max(nSnp),
#     inbred = TRUE,
#     species = "WHEAT"
#   )
# }

founderPop = runMacs2(
  nInd = nParents,
  nChr = nChr,
  inbred = TRUE,
  mutRate=2.5e-10,
  ploidy=2L
)
print(founderPop)

SP = SimParam$new(founderPop)

# Add SNP chip/array (with markers that don't overlap with QTL)
SP$restrSegSites(nQtl*4, max(nSnp))
if (any(nSnp > 0)) {
  for (n in nSnp) {
    SP$addSnpChip(n)
  }
}

# Add a trait such as yield (effects sampled from a Gaussian distribution)
# Trait 1: Additive, and GxE, gaussian
SP$addTraitAG(
  nQtlPerChr = nQtl,
  mean = initMeanG,
  var = initVarG,
  varEnv = initVarEnv,
  varGxE = initVarGE
)


# Trait 2: Additive, epistatic, and GxE, gaussian
SP$addTraitAEG(
  nQtlPerChr = nQtl,
  mean = initMeanG,
  var = initVarG,
  relAA = initRelAA,
  varEnv = initVarEnv,
  varGxE = initVarGE
)


# Add another trait (effects sampled from a gamma distribution)
# Trait 3: Additive, and GxE, gamma
SP$addTraitAG(
  nQtlPerChr = nQtl,
  mean = initMeanG,
  var = initVarG,
  varEnv = initVarEnv,
  varGxE = initVarGE,
  gamma = TRUE,
  shape = 0.25
)
# Trait 4: Additive, epistasis, and GxE, gamma
SP$addTraitAEG(
  nQtlPerChr = nQtl,
  mean = initMeanG,
  var = initVarG,
  relAA = initRelAA,
  varEnv = initVarEnv,
  varGxE = initVarGE,
  gamma = TRUE,
  shape = 0.25
)

source(file = paste(here(), breedingProgDir, "ExportQtlEffects.R", sep="/"))

# Collect pedigree
SP$setTrackPed(TRUE)

# Create founder parents
Parents = newPop(founderPop)

# Add phenotype reflecting evaluation in EYT
Parents = setPheno(Parents, varE = varE, reps = repEYT)

rm(founderPop)
