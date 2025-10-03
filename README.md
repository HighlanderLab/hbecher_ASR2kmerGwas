# hbecher_ASR2kmerGwas
Code to convert ASR simulations to k-mer datasets

## Workflow
1. Simulate sequences with ASR (`01haplogenomesExample.R`)
2. Convert to FASTA files (`01haplogenomesExample.R`)
3. Generate sequencing reads from the FASTA files (e.g. wgsim, `02makeReads.sh`)
4. Count k-mers (e.g. with kmc, `03countKmers.sh`)

## Requirements
- R/AlphaSimR
- wgsim (https://github.com/lh3/wgsim/)
- kmc (https://anaconda.org/bioconda/kmc)