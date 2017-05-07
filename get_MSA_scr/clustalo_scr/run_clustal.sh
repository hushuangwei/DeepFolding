Seq_Path="/home/xubiao/scratch/blast_test/extract_seq/mview-master/seqs";
output="clustalo_output"
mkdir -p $output
readarray arr < $1
for i in ${arr[@]}
do
name=$(echo $i | tr '[:upper:]' '[:lower:]')
name1=${Seq_Path}/$name;
if [ -e $name1 ];then
echo $name1
./clustalo-1.2.4-Ubuntu-x86_64 -i $name1 -o ${output}/$name
fi
done