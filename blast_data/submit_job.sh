name=$(head -500 PDB_list.txt)
for i in $name
do
name1=$(echo $i | tr '[:upper:]' '[:lower:]')
echo $name1
cd result_1_500
qsub -v pdbname=$name1 ../run_gromacs2.pbs
cd ..
done