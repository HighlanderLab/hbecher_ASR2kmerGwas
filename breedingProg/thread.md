As promised I took scripts of our past  wheat simulation  with AlphaSimR (see https://github.com/HighlanderLab/jbancic_alphasimr_plants/blob/main/01_LineBreeding/LineProgram.jpeg) and adapted it to the needs of our collaboration.

The scripts are at https://www.dropbox.com/scl/fo/d3i2zvs2v8zhczd0brd1e/ACgnDIFSaI2LDCEn-t4ron0?rlkey=r3gem5i6o7zp0cw6ig51klcix&dl=0. This is all R script code and AlphaSimR can be installed from CRAN. It helps if your R installation is OpenMP enabled so we can parallelise some steps in the code (done automatically).

I am happy to walk you through the code in the meeting, but very quickly some points:

    00RUNME.R script runs the simulation, including export of the data for your analysis
    GlobalParameters.R script defines various simulation parameters, such as:
        number of years a breeding programme is run (will impact the size of the data we generate)
        how many markers we capture
        how many loci affect a trait
        proportion of variation due to additive and epistasis effects
        how many individuals we generate at each stage (see the above scheme LineProgram.jpeg)
        …

I assume you will want to run this simulation in your environment so we are not moving output files around. These scripts will generate these files for your analysis:

    snpGenoChipXTrain.txt is holding id column + nInd * nSnp genotypes for training individuals for snp array/chip X (defined by nSnp variable in GlobalParameters.R)
    phenoTrain.txt is holding id column and phenotype value column for training individuals
    gvTrain.txt is holding id column, genetic value column, breeding value column (only statistical additive component of genetic value), and epistatic value column (only statistical epistatic component of genetic value) for training individuals
    snpGenoChipXPredictY.txt is holding id column + nInd * nSnp genotypes for prediction individuals for snp array/chip X (defined by nSnp variable in GlobalParameters.R) - we have several stages in a breeding programme where we might want to do prediction hence we have different Y populations
    phenoPredictY.txt is holding id column and phenotype value column for prediction individuals in population Y (we don't really need this, but got to have just in case)
    gvPredictHDRW.txt is holding id column, genetic value column, breeding value column (only statistical additive component of genetic value), and epistatic value column (only statistical epistatic component of genetic value) prediction individuals in population Y - we want to compare our model predictions against these values

I can also show how to get access to simulated additive and epistatic effects for causal loci - see CreateGenomeTraitsParents.R script how to access these quantities.



------------------------------------------------------
I have now updated my simulation at https://www.dropbox.com/scl/fo/d3i2zvs2v8zhczd0brd1e/ACgnDIFSaI2LDCEn-t4ron0?rlkey=r3gem5i6o7zp0cw6ig51klcix&dl=0 

Specifically:

    I added ExportQtlEffects.R to export information which of the tracked SNP are causal (are QTL=Quantitative Trait Loci), what are their additive effects, which other QTL they interact with, and what is the effect of that epistatic interaction
    I tweaked the GlobalParameters.R to run burn-in of the breeding programme for 10 years (generates additional structure in genome due to selection), but I only collect training data in year 10 (this will give us a bit more than 1000 accessions, which you can further subset if you want, or decrease the nPYT variable in the script)
    I tweaked CreateGenomeTraitsParents.R such that we now have one trait that has QTL effects sampled from a Gaussian distribution and another trait that has QTL effects sampled from a Gamma distribution - will show the impact below

Files of interest for QTL (this is ground truth - we don't know this in real life):

    QTL1Gauss.txt - QTL, chromosome, locus, additive effect, QTL2 (in interaction), chromosome2, locus2, epistatic effect of interaction - all for trait1 (Gauss QTL effects)
    QTL2Gamma.txt - as above but for trait1 (Gamma QTL effects)
    QTL1GaussAddEff.png, QTL1GaussEpiEff.png, QTL2GammaAddEff.png, QTL2GammaEpiEff.png - distributions of QTL effects
    QTL1GaussAddEffMap.png, QTL1GaussEpiEffMap.png, QTL2GammaAddEffMap.png, QTL2GammaEpiEffMap.png - map of QTL effects along the genome

I attach examples. Note that the "epistatic map" QTL*EpiEffMap.png can be misleading - it's impossible to map pairwise epistatic interactions onto one dimension, since we have two QTL that interact and they can be close or far apart from each other. I created a 2D version as well and tried to denote size of epistatic effect by the size of the dot (not easy!) - see the attached.

Files of interest for individuals (for modelling):

    phenoTrain.txt - individual, phenotypic value for trait1 (Gauss QTL effects), and phenotypic value for trait2 (Gamma QTL effects) - this is for training individuals
    gvTrain.txt - individual, genetic value (gv) for trait1 and trait2, breeding value (bv) for trait1 and trait2, and epistatic deviations (aa) for trait1 and trait2 - this is ground truth upon which phenoTrain.txt is generated (we don't know this in real life)
    snpGenoChip1Train.txt - individual and SNP genotypes for SNP chip/array 1
    phenoPredictHDRW.txt - as above but for a new generation (set of individuals called HDRW, among others), which we would like to predict
    gvPredictHDRW.txt - as above but for a new generation, which we would like to predict
    snpGenoChip1PredictHDRW.txt - as above but for a new generation, which we would like to predict

As you can see in the attached, the simulation generates extremely complex polygenic traits - each impacted by 16,000 QTL and 8,000 pairwise interactions! This is a very very hard case for GWAS, particularly if we only have ~200 accessions. We have options here - further decrease the shape parameter in the gamma distribution to get heavier tails and/or reduce the number of QTL!

Let me know if you can run the code and if you have any questions.
