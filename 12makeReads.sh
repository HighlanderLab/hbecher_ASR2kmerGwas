

wgsim="wgsim/wgsim"  # path to wgsim binary (code obtained from https://github.com/lh3/wgsim and compiled)



inDir="haploGenomes"


rl=150 # read length
dep=50 # sequencing depth (per haploid genome size)
gs=10000 # genome size (haploid)
for f in `ls $inDir/ind*.fasta`; do

    base=$(basename $f .fasta)
    # -e is the base error rate
    # -d is the avg. outer distance between the two ends
    # -s is the stddev of the outer distance
    # -1 and -2 are the read lengths
    # -N is the number of read pairs to simulate
    echo "Simulating reads for $base"
    $wgsim -r 0 -R 0 -e 0.0001 -d 550 -s 50  \
          -1 $rl -2 $rl -N $((gs/2/$rl*$dep)) \
          $f haploGenomes/${base}_fw.fq  haploGenomes/${base}_rw.fq    > /dev/null
done