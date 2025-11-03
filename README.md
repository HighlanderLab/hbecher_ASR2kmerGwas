# hbecher_ASR2kmerGwas
Code to convert ASR simulations to k-mer datasets

## Workflow
1. Run `Rscript 10ASR2FASTA.R`.
2. This will produce
  - 3000 `FASTA` files (multi-chromsome, diploid genomes)
  - `AllGvData.txt` (tab-separated, 1st column matches the FASTA file names)
  - `AllPhenoData.txt` (tab-separated)

## For k-mer GWAS
Use column `gv.Trait1` from file `AllGvData.txt` as phenotypes.

## Requirements
- R: AlphaSimR, dplyr, here
