 clc
 clear all
% list=importdata('list');
aa = containers.Map;
aa('GLY') = 'G';
aa('ALA') = 'A';
aa('ILE') = 'I';
aa('LEU') = 'L';
aa('PHE')='F';
aa('PRO')='P';
aa('GLN')='Q';
aa('GLU')='E';
aa('ARG')='R';
aa('ASN')='N';
aa('ASP')='D';
aa('CYS')='C';
aa('HIS')='H';
aa('LYS')='K';
aa('MET')='M';
aa('SER')='S';
aa('THR')='T';
aa('TRP')='W';
aa('TYR')='Y';
aa('VAL')='V';
PDB_struct = pdbread('PDBs_new/1bg6_1_A.pdb');
PDB_struct_ref = pdbread('PDBs/1bg6.pdb');
seq_ref=PDB_struct_ref.Sequence.Sequence;
count_seq=0;
for i=1:12
    if strcmp(PDB_struct.Model.Atom(i).AtomName,'CA') 
        count_seq=count_seq+1;
    resname=PDB_struct.Model.Atom(i).resName
    coor_seq_init(count_seq)=aa(resname);
    end
end
k = strfind(seq_ref,coor_seq_init)