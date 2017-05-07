#name=$(awk '{if (NR>500 && NR<=1000) print $0}' PDB_list.txt)
#for i in $name
#do
#name1=$(echo $i | tr '[:upper:]' '[:lower:]')
#echo $name1
#cd result_1_500
#qsub -v pdbname=$name1 ../run_gromacs2.pbs
#cd ..
#done
readarray name_array < $1
length=${#name_array[@]}
echo $length
for ((i=0;i<$length;i=i+100))
do
echo $i
#cd result_3000_6000 
#cd result_rest 
sub_name=${name_array[@]:$i:100}
#echo $i
#echo $sub_name
pdbname="${sub_name[@]}";
#name1="$(echo $pdbname | tr '[:upper:]' '[:lower:]')" 
#echo ${name1[1]}
qsub -v pdbname="${sub_name[@]}" run_mview.pbs
#cd ..
#echo $name1
#echo "-------------------"

done
