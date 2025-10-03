

wgsim="wgsim/wgsim"  # path to wgsim binary (code obtained from https://github.com/lh3/wgsim and compiled)



inDir="haploGenomes"


rl=150 # read length
dep=50 # sequencing depth (per haploid genome size)
gs=10000 # genome size (haploid)
for f in `ls $inDir/ind*.fasta`; do

    base=$(basename $f .fasta)
    $wgsim -r 0 -R 0 -e 0.0001 -d 550 -s 50  \
          -1 $rl -2 $rl -N $((gs/2/$rl*$dep)) \
          $f haploGenomes/${base}_fw.fq  haploGenomes/${base}_rw.fq    > /dev/null
done