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
  mutRate=2e-10,
  ploidy=2L
)
print(founderPop)

SP = SimParam$new(founderPop)

# Add SNP chip/array (with markers that don't overlap with QTL)
SP$restrSegSites(nQtl, max(nSnp))
if (any(nSnp > 0)) {
  for (n in nSnp) {
    SP$addSnpChip(n)
  }
}

# Add a trait such as yield (effects sampled from a Gaussian distribution)
SP$addTraitAEG(
  nQtlPerChr = nQtl,
  mean = initMeanG,
  var = initVarG,
  relAA = initRelAA,
  varEnv = initVarEnv,
  varGxE = initVarGE
)

# Add another trait (effects sampled from a gamma distribution)
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
