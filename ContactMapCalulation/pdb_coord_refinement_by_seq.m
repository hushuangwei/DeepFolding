clc
clear all
% PDB.SrcPath = 'PDBs';
% PDB.SrcFile = dir(fullfile(PDB.SrcPath, '*.pdb'));
% PDB.SrcFile = {PDB.SrcFile.name};
% % %
% % %
% list=PDB.SrcFile;
!mkdir -p coords_refined
%problematic_pdbs=[];
list=importdata('total_list_test');
count=0;
count_good=0;
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



% %

for pdbid=1:length(list)
PDB_struct = pdbread(list{pdbid});
pdbname=list{pdbid}
pdbname=pdbname(10:13)
PDB_struct_ref = pdbread(['PDBs/' pdbname '.pdb']);
seq_ref=PDB_struct_ref.Sequence.Sequence;
coor_seq_init=[];
if length(PDB_struct.Model.Atom)<8
    count=count+1;
    problematic_pdbs{count}={pdbname};
    continue;
end    
count_seq=0;
skip_resnum=0;
i=1;
while i<=8
    if strcmp(PDB_struct.Model.Atom(i+skip_resnum).AtomName,'CA') 
        count_seq=count_seq+1
        resname=PDB_struct.Model.Atom(i+skip_resnum).resName
    if  ~isKey(aa,resname) || PDB_struct.Model.Atom(i+skip_resnum).resSeq-PDB_struct.Model.Atom(i+skip_resnum).resSeq~=0
        skip_resnum=skip_resnum+1;
        count_seq=0;
        i=1;
        continue;
    end
%     resname
%     aa(resname)
    coor_seq_init(count_seq)=aa(resname)
    end
    i=i+1;
end
seq_ref
coor_seq_init
k_pos = strfind(seq_ref,coor_seq_init)
if (isempty(k_pos))
    count=count+1;
    problematic_pdbs{count}={pdbname};
    continue;
end

resid_start=PDB_struct.Model.Atom(1).resSeq-k_pos(1)+1;
break_flag=0;
%resid_start=PDB_struct.Model.Atom(1).resSeq;
resid_end=PDB_struct.Model.Atom(end).resSeq;
if k_pos(1)==1 && resid_end-resid_start+1==length(seq_ref)
    pdbname
    count_good=count_good+1;
    good_pdbs{count_good}={pdbname};
end
% if resid_end<resid_start
%     count=count+1;
%     problematic_pdbs{count}={pdbname};
%     continue;
% end
%coord=nan(resid_end-resid_start+1,3);
coord=nan(length(seq_ref),3);
for i=1:length(PDB_struct.Model.Atom)
   if strcmp(PDB_struct.Model.Atom(i).AtomName,'CA')
       resid=PDB_struct.Model.Atom(i).resSeq;
       if resid<resid_start
       break_flag=1;   
       break;
       end
   coord(resid-resid_start+1,1)=PDB_struct.Model.Atom(i).X;
   coord(resid-resid_start+1,2)=PDB_struct.Model.Atom(i).Y;
   coord(resid-resid_start+1,3)=PDB_struct.Model.Atom(i).Z;
   end
end
if break_flag==1
    count=count+1;
    problematic_pdbs{count}={pdbname};
    continue;
end
% D=CalculateDistanceMatrix(coord);
 if any(isnan(coord))
    dlmwrite(['coords_refined/' pdbname '_star.txt'],coord,'delimiter','\t')
 else
    dlmwrite(['coords_refined/' pdbname '.txt'],coord,'delimiter','\t')
 end
%break
end