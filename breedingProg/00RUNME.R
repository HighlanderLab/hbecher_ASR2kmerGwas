# setwd(dir = "~/Storages/DropBox/5_Transfer/WheatSimEpistasis")

# Phenotypic wheat breeding program with doubled haploid technology and
# collecting genomic data

# ---- Clean environment and load packages ----
# rm(list = ls()) # moved to main script
# install.packages(pkgs = "AlphaSimR")
#library(package = "AlphaSimR")

# ---- Load global parameters ----
source(file = here(breedingProgDir, "GlobalParameters.R"))
scenarioName = "LinePheno_DH"

# ---- Create list to store results from reps ----
results = list()

for (REP in 1:nReps) {
  # REP = 1
  cat("Working on REP:", REP, "\n")

  # ---- Create a data frame to track key parameters ----
  output = data.frame(
    year = 1:nCycles,
    rep = rep(REP, nCycles),
    scenario = rep(scenarioName, nCycles),
    meanG = numeric(nCycles),
    varG = numeric(nCycles),
    accSel = numeric(nCycles)
  )

  # ---- Create initial parents ----
  source(file = here(breedingProgDir, "CreateGenomeTraitsParents.R"))

  # ---- Fill breeding pipeline with unique individuals from initial parents ----
  source(file = here(breedingProgDir, "FillPipeline.R"))

  # ---- Burn-in phase ----
  for (year in 1:nBurnin) {
    # year = 1
    cat("  Working on burn-in year:", year, "\n")
    source(file = here(breedingProgDir, "UpdateParents.R")) # Pick parents
    source(file = here(breedingProgDir, "AdvanceYear.R")) # Advances yield trials by a year
    source(file = here(breedingProgDir, "StoreTrainPop.R")) # Store training population
    # Report results
    output$meanG[year] = meanG(DH)[1]
    output$varG[year] = varG(DH)[1, 1]
  }

  # Save genotype data for external use
  # for (chip in 1:length(nSnp)) {
  #   # chip = 1
  #   snpGeno = pullSnpGeno(TrainPop, snpChip = chip)
  #   write.table(
  #     x = snpGeno,
  #     file = paste0("snpGenoChip", chip, "Train.txt"),
  #     sep = " ",
  #     col.names = FALSE, # TODO: uncomment this if you want column names
  #     row.names = TRUE,
  #     quote = FALSE
  #   )
  # }

  # Save phenotype data for external use
  phenos = data.frame(id = TrainPop@id, pheno = pheno(TrainPop))
  # head(phenos)
  write.table(
    x = phenos,
    file = paste0("phenoTrain", year, ".txt"),
    sep = " ",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )

  # Save genetic values for external use
  tmp = data.frame(
    id = TrainPop@id,
    gv = gv(TrainPop),
    bv = bv(TrainPop),
    aa = aa(TrainPop)
  )
  # head(tmp)
  write.table(
    x = tmp,
    file = paste0("gvTrain", year, ".txt"),
    sep = " ",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )

  # ---- Future phase ----
  for (year in (nBurnin + 1):(nBurnin + nFuture)) {
    # year = nBurnin + 1
    cat("  Working on future year:", year, "\n")
    source(file = here(breedingProgDir, "UpdateParents.R")) # Pick parents
    source(file = here(breedingProgDir, "AdvanceYear.R")) # Advances yield trials by a year
    # Report results
    output$meanG[year] = meanG(DH)[1]
    output$varG[year] = varG(DH)[1, 1]

    if("TrainPop" %in% ls()) rm(TrainPop) # Remove training population object to save memory
    for (pop in c("HDRW", "PYT", "AYT", "EYT")) {
      # pop = "HDRW"
      # Save genotype data for external use
      # for (chip in 1:length(nSnp)) {
      #   # chip = 1
      #   snpGeno = pullSnpGeno(get(pop), snpChip = chip)
      #   write.table(
      #     x = snpGeno,
      #     file = paste0("snpGenoChip", chip, "Predict", pop, ".txt"),
      #     sep = " ",
      #     col.names = FALSE, # TODO: uncomment this if you want column names
      #     row.names = TRUE,
      #     quote = FALSE
      #   )
      # }

      haplots = pullSegSiteHaplo(get(pop), simParam = SP)
      ff <- gzfile(paste0("haplotypes", pop, "yesr", year, ".txt.gz"), "wt")
        write.table(
          x = haplots,
          file = ff,
          sep = " ",
          col.names = FALSE, # TODO: uncomment this if you want column names
          row.names = TRUE,
          quote = FALSE
        )
        close(ff)
      
      # Save phenotype data for external use
      phenos = data.frame(id = get(pop)@id, pheno = pheno(get(pop)))
      # head(phenos)
      write.table(
        x = phenos,
        file = paste0("phenoPredict", pop, "year", year, ".txt"),
        sep = " ",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )

      # Save genetic values for external use
      tmp = data.frame(
        id = get(pop)@id,
        gv = gv(get(pop)),
        bv = bv(get(pop)),
        aa = aa(get(pop))
      )
      # head(tmp)
      write.table(
        x = tmp,
        file = paste0("gvPredict", pop, "year", year, ".txt"),
        sep = " ",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }
  }

  # Save results from current replicate
  results = append(results, list(output))
}

# Save results
#saveRDS(results, file = paste0(scenarioName, ".rds"))

# ---- Analyze results ----
#source(file = here(breedingProgDir, "ANALYZERESULTS.R"))
