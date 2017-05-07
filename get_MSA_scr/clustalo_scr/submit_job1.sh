readarray name_array < $1

for ((i=0;i<100;i=i+100))
do
echo $i
sub_name=${name_array[@]:$i:100}
#echo $sub_name
pdbname="${sub_name[@]}";
name1="$(echo $pdbname | tr '[:upper:]' '[:lower:]')" 
echo $name1
#echo ${#name1[@]}
qsub  -v pdbname="$(echo ${sub_name[@]} | tr '[:upper:]' '[:lower:]')" run_clustal.pbs
done