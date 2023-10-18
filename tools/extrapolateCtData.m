function [sinogram, CtData] = extrapolateCtData(D, allAngles)
%EXTRAPOLATECTDATA Extrapolates real data from CtData structure to given
%random angles and creates a new CtData structure
%
%   INPUT
% D         Original data with fixed set of projection angles
% allAngles T sequences of random projection angles used for the simulation
%
%   OUTPUT
% sinogram  Sequence of sinograms WITHOUT any added noise
% CtData    HelTomo structure for creating the (block diagonal) forward operator
%
% T H   2023

T = length(D);

Nang = size(allAngles,2);
switch size(allAngles,1)
    case 1
        % Use the same angles for all time steps
        allAngles = repmat(allAngles,T,1);
    case T
        % Do nothing
    otherwise
        error('Incorrect number of angle-vectors given: %d when %d are needed!',size(allAngles,1),T);
end

CtData = D;
Ndet = D(1).parameters.numDetectors;

sinogram = zeros(Nang,Ndet,T);

if max(abs(allAngles(:))) > 360
    warning('Maximum angle is greater than 360! This could be problematic.')
end

for t = 1:T
    CtData(t).sinogram = []; % No need for this
    OGsino = [D(t).sinogram; D(t).sinogram(1,:)];
    OGangles = [D(t).parameters.angles, D(t).parameters.angles(1) + 360];

    angles = allAngles(t,:);
    CtData(t).parameters.angles = angles;
    CtData(t).numberImages = Nang;

    % Integer (shifted by 1) and fractional part of the random angles
    angInd = floor(angles) + 1;
    angWght = rem(angles(:),1); % Forced to column vector!

    if ~all(ismember(angInd,OGangles)) || ~all(ismember(angInd-1,OGangles))
        error('Can not extrapolate because projection is missing!')
    end

    sinogram(:,:,t) = angWght.*OGsino(angInd,:) + (1-angWght).*OGsino(angInd+1,:);
end