# hbecher_ASR2kmerGwas

Code to convert ASR simulations to k-mer datasets

## Workflow

1.  Run `Rscript 10ASR2FASTA.R`.
2.  This will produce

-   3000 `FASTA` files (multi-chromsome, diploid genomes)
    -   For more files, up to 33000, comment out line 126 (`haploFiles <- dir(pattern="haplotypes*PYT")`) of file `10ARS2FASTA.R`)
-   `AllGvData.txt` (tab-separated, 1st column matches the FASTA file names)
-   `AllPhenoData.txt` (tab-separated)
-   `AllVarLociAddPart.txt` (tab-separated) with one column for each trait, showing the additive effect tat each variant site ahs on each trait
-   Four per-trait QTL files. These also list the pairwise epistatic interaction effect for each trait.

## For k-mer GWAS

Use column `gv.Trait1` from file `AllGvData.txt` as phenotypes.

## Requirements

-   R: AlphaSimR, dplyr, here

## Code sources

This uses code from <https://github.com/Rabab53/ggenerator>.
