function D = CalculateDistanceMatrix(Coords, I, J, SwitchBSX)
%% CALCULATEDISTANCEMATRIX 
% Description: calculate the distance matrix of atoms
% Author: Jin Dai
% Created Date: 2014.06.04
% Last Modified Date: 2016.09.26


%% Default Input Settings
% 1 input argument {Coords} : return a matrix of the distance between all atoms

% 2 input argument {Coords, I} : return a vector of the distance from the atom specified by I to all atoms

% 3 input argument {Coords, I, J} : return a scalar of the distance from the atom specified by I to the atom specified by J

% 4 input argument {Coords, I, J, SwitchBSX} : switch to control the code version calculation

if nargin < 4
    SwitchBSX = true;
end

if nargin < 3 || isempty(J)
    J = 'ALL';
end

if nargin < 2 || isempty(I)
    I = 'ALL';
end

if nargin < 1
    fprintf('Please input one coordinate-vector as the argument!\n');
    return;
end

[Row Column] = size(Coords);
if Column ~= 3
    Coords = Coords';
    [Row Column] = size(Coords);
end


%% Calculate by matrix method
Mij = Coords*Coords';% get coupled matrix Mij = Coords*Coords'
DiagVec = diag(Mij); % diagonal element of Mij to form a column vector

if SwitchBSX
    DiagMtr = DiagVec*ones(1, length(DiagVec)); % expand this vector into a matrix as {DiagVec(1) DiagVec(2) ... }
    % get DistanceMatrix = sqrt(DiagMtr + DiagMtr' - 2*Mij), DistanceMatrix is the diastance matrix between all atoms
    DistanceMatrix = sqrt(DiagMtr + DiagMtr' - 2*Mij);
else
    % get DistanceMatrix by bxsfun
    DistanceMatrix = sqrt(bsxfun(@plus, DiagVec, DiagVec') - 2*Mij);
end


%% return D depending on the input arguments
if strcmp(J, 'ALL')
    if strcmp(I, 'ALL')
        D = DistanceMatrix;
    else
        D = DistanceMatrix(:, I);
    end
else
    D = DistanceMatrix(J, I);
end


end

