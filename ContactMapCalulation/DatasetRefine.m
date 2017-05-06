clear all;
clc;

RefreshFiles = false;


%%
PDB.SrcPath = 'PDBs';
PDB.SrcFile = dir(fullfile(PDB.SrcPath, '*.pdb'));
PDB.SrcFile = {PDB.SrcFile.name};

PDB.ObjPath = [PDB.SrcPath '_Refine'];
if isempty(dir(PDB.ObjPath))
    mkdir(PDB.ObjPath);
elseif (RefreshFiles == true)
    rmdir(PDB.ObjPath, 's');
    mkdir(PDB.ObjPath);
end
PDB.ObjFile = cell(size(PDB.SrcFile));

PDB.ObjCell = cell(size(PDB.SrcFile));
PDB.ObjList = fullfile('.', ['RefineList_' PDB.SrcPath '.txt']);
%%
list=importdata('list');

%%
% print information and start timer
fprintf('Refining pdb dataset starts...\n');
tic;

% temporary variables
TempSrcPath = PDB.SrcPath;
TempSrcFile = PDB.SrcFile;
TempObjPath = PDB.ObjPath;
TempObjCell = PDB.ObjCell;
% %
% start parpool
TempPool = gcp('nocreate');
if isempty(TempPool)
    switch computer('arch')
        case {'win64', 'win32'}
            TempPoolNum = 2;
        case 'glnxa64'
            TempPoolNum = 4;
        case 'maci64'
            TempPoolNum = 1;
        otherwise
            TempPoolNum = 1;
    end
    TempPool = parpool('local', TempPoolNum);
end
%
% main loop
for iF = 1:numel(list)
%     if iF<=6143
%         continue
%     end
list{iF}
    TempStruct = pdbreadMod(fullfile(TempSrcPath, list{iF}));
    
    TempObjCell{iF} = sprintf('%s :', list{iF});
    %--
    %fprintf(TempFidList, '%s :', TempSrcFile);
    
    for iM = 1:numel(TempStruct.Model)
        % check model
        TempModel = TempStruct.Model(iM);
        if ~isfield(TempModel, 'Atom')
            TempObjCell{iF} = [TempObjCell{iF}  sprintf(' %s', 'NULL')];
            %--
            %fprintf(TempFidList, ' %s', 'NULL');
            continue;
        end
        
        % check chain
        TempChainId = unique({TempModel.Atom.chainID});
        for iC = 1:numel(TempChainId)
            % file name
            TempSrc = list{iF};
            TempFileName = [TempSrc(1:4) '_' num2str(iM) '_' TempChainId{iC} '.pdb'];
            % print information in list
            TempObjCell{iF} = [TempObjCell{iF} sprintf(' %s', TempFileName)];
            %--
            %fprintf(TempFidList, ' %s', TempFileName);
            
            % get chain index
            TempIdxChain = strcmp({TempModel.Atom.chainID}, TempChainId{iC});
            
            % get CA index
            TempIdxCA = strcmp({TempModel.Atom.AtomName}, 'CA');
            % get CB index
            TempIdxCB = strcmp({TempModel.Atom.AtomName}, 'CB');
            % get HA3 index
            TempIdxHA3 = strcmp({TempModel.Atom.AtomName}, 'HA3');
            
            % index for AtomName
            TempIdxAtomName = find((TempIdxCA | TempIdxCB | TempIdxHA3) & TempIdxChain);
            
            % atoms in specific chain
            TempAtom = TempModel.Atom(TempIdxAtomName);
            
            % -----------------------------------------------------------------
            % TempModel -> TempAtom : array index changes
            % -----------------------------------------------------------------
            % shift resSeq by iCode [loop-version]
            TempResAll = strcat(cellfun(@num2str, {TempAtom.resSeq}, 'UniformOutput', false), {TempAtom.resName}, {TempAtom.iCode});
            TempResList = unique(TempResAll, 'stable');
            
            TempShift = 0; % shift number of residue number
            TempKeep = zeros(numel(TempAtom),1); % marker for keepable atoms
            for iR = 1:numel(TempResList)
                % atom index in residue
                TempIdxAtomInRes = find(strcmp(TempResAll, TempResList{iR}));
                % increase shift if iCode is not empty
                if ~strcmp(TempAtom(TempIdxAtomInRes(1)).iCode, '')
                    TempShift = TempShift + 1;
                end
                
                % atom name in residue
                TempAtomNameInRes = {TempAtom(TempIdxAtomInRes).AtomName};
                
                % check altLoc for CA
                TempJ = find(strcmp(TempAtomNameInRes, 'CA'));
                if ~isempty(TempJ)
                    TempOccupancy = [TempAtom(TempIdxAtomInRes(TempJ)).occupancy];
                    TempMaxOcp = TempJ(TempOccupancy == max(TempOccupancy));
                    TempKeep(TempIdxAtomInRes(TempMaxOcp(1))) = 1;
                end
                
                % check altLoc for CB or HA3
                TempJ = find(strcmp(TempAtomNameInRes, 'CB') | strcmp(TempAtomNameInRes, 'HA3'));
                if ~isempty(TempJ)
                    TempOccupancy = [TempAtom(TempIdxAtomInRes(TempJ)).occupancy];
                    TempMaxOcp = TempJ(TempOccupancy == max(TempOccupancy));
                    TempKeep(TempIdxAtomInRes(TempMaxOcp(1))) = 1;
                end
                
                % reassign the resSeq and iCode for insertion residues
                for iA = 1:numel(TempIdxAtomInRes)
                    TempI = TempIdxAtomInRes(iA);
                    TempAtom(TempI).resSeq = TempAtom(TempI).resSeq + TempShift;
                    if ~strcmp(TempAtom(TempI).iCode,'')
                        TempAtom(TempI).iCode = 'X';
                    end
                end
                
            end
            
            % shift resSeq by iCode [vector-version]
            % ??? How ???
            
            % record results
            Target = struct('Model', []);
            Target.Model.Atom = TempAtom(logical(TempKeep));
            %Target.Model.Terminal = CreateEmptyTerminal();
            
            % write to pdb file
            warning('off', 'bioinfo:pdbwrite:MissingAllMandatoryFields');
            pdbwrite(fullfile(TempObjPath, TempFileName), Target);
        end
    end
end

% stop parpool
if ~isempty(TempPool)
    delete(TempPool);
end

% record result in PDB.ObjCell
PDB.ObjCell = TempObjCell;

% print and stop timer
fprintf('Refining pdb dataset completed!\n');
toc;
fprintf('\n');


%% 
TempFidList = fopen(PDB.ObjList, 'w');
if TempFidList == -1
    fprintf('ERROR: unable to open file ''%s''!\n', PDB.ObjList);
    return;
end

for iO = 1:numel(PDB.ObjCell)
    fprintf(TempFidList, '%s\n', PDB.ObjCell{iO});
end

fclose(TempFidList);


%%
clear Temp*;
clear i*;
clear ans;

