for f in $1/psiblastDsbAOut*
do
#name1=${f%*/}
#name=${name1#*_}
name=${f:(-4)}
echo $name
mview $f -out fasta | sed 's/-//g' > seqs/$name
done