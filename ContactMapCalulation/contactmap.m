clear
clc
!mkdir -p contactmap_result_CA_8A
PDBs = dir(fullfile('result_test','*.txt'));
PDBname = {PDBs.name};
% %
for iF = 1:numel(PDBname)
A=dlmread(['result_test/' PDBname{iF}]);
B=(A<8);
dlmwrite(['contactmap_result_CA_8A/' PDBname{iF}],B,'delimiter','\t')
end