---
editor_options: 
  markdown: 
    wrap: 72
---

# hbecher_ASR2kmerGwas

Code to convert ASR simulations to k-mer datasets

## Workflow

1.  Run `Rscript 10ASR2FASTA.R`.
2.  This will produce:

-   One reference genome FASTA file `fatsa/ReferenceGenome.fa`. All
    nucleotide positions in other files are given relative to this file.
-   3000 `FASTA` files of simulated individuals (multi-chromsome,
    diploid genomes)
    -   For more files, up to 33000, comment out line 126
        (`haploFiles <- dir(pattern="haplotypes*PYT")`) of file
        `10ARS2FASTA.R`)
-   `AllGvData.txt` (tab-separated, 1st column matches the FASTA file
    names)
-   `AllPhenoData.txt` (tab-separated)
-   `AllVarLociAddPart.txt` (tab-separated) with one column for each
    trait, showing the additive effect that each variant site has on
    each trait. This file can also be used to match variant IDs (1st
    column) to nucleotide positions in the reference (cols 2 and 3).
-   Four per-trait QTL files. These also list the pairwise epistatic
    interaction effect for each trait.

## Traits

Four traits ares simulated:

-   Trait1: purely additive, gaussian
-   Trait2: additive + epistatic, gaussian
-   Trait3: purely additive, gamma
-   Trait4: additive + epistatic, gamma

The additive contributions of each variant site (position with respect
to the reference genome), are written to `AllVarLociAddPart.txt`. For
the epistatic traits, the pairwise interaction effects are written to
the per-trait QTL files, `QTL[...]txt`.

## For k-mer GWAS

Use one of the `gv.Trait...` columns from file `AllGvData.txt` as the
phenotypes. The do not have any environmental randomness. For more
environmental noise, you can use `pheno.Trait...` from file
`AllPhenoData.txt` instead.

## Requirements

-   R: AlphaSimR, dplyr, here

## Code sources

This uses code from <https://github.com/Rabab53/ggenerator>.
