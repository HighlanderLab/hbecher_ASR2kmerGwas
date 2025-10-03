
# make sure the k-mer toolkit kmc is in path
# can be installed via anacond
source /home/hbecher/miniconda3/etc/profile.d/conda.sh
conda activate kmc

inDir="haploGenomes"
for i in `ls $inDir/ind*.fasta`; do
    echo "$i"
    base=$(basename $i .fasta)
    echo $base
    # make file listing the input read files for each individual
    printf "$inDir/$base"_fw".fq\n$inDir/$base"_rw".fq" > $inDir"/readfiles_"$base
    echo "Counting k-mers for $base"
    kmc -k21 -t3 -cs500000000 -fq -ci1 -fq @$inDir/readfiles_$base $inDir/${base}_db21 $inDir
    echo "Generating spectrum for $base"
    kmc_tools transform $inDir/${base}_db21 histogram $inDir/${base}_hist21.txt
    #rm $inDir/readfiles_$base
done

