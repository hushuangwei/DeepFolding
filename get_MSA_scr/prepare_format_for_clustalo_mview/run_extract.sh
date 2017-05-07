readarray name_array < $1
length=${#name_array[@]}
echo $length
sub_name=${name_array[@]:300:100}
pdbname="${sub_name[@]}";
name1=($pdbname)
for i in {0..99}
do
f=${name1[$i]}
echo $f
name=${f:(-4)}
echo $name
mview $f -out fasta | sed 's/-//g' > seqs/$name
done


