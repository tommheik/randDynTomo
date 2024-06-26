function [sinogram, CtData] = extrapolateStempoCtData(D, allAngles)
% EXTRAPOLATESTEMPOCTDATA Extrapolates stempo data from CtData structure to
% given random angles and creates a new CtData structure
%
%   INPUT
% D         Original data with fixed set of projection angles
% allAngles T sequences of random projection angles used for the simulation
%           The angles should go beyond 360 degrees for t > 2 at this point
%
%   OUTPUT
% sinogram  Sequence of sinograms WITHOUT any added noise
% CtData    HelTomo structures for creating the (block diagonal) forward operator
%
% T. Heikkilä   2023

% Angles should be generated by following formula:
% shifts = linspace(0,8*360 - 188,T); % 8*360 - 8 is maximum angle and we
%                                       sample from [0, 180] degrees
% allAngles = 180*rand(1,23) + shifts' % 23 random angles from s + [0, 180]
% OR have same shift but for last time step only sample to 172 deg!


T = size(allAngles,1);
Nang = size(allAngles,2);

CtData = D;
OGsino = D.sinogram;
Ndet = size(OGsino,2);

sinogram = zeros(Nang,Ndet,T);
angReso = D.parameters.angles(2) - D.parameters.angles(1); % Angular resolution

for t = 1:T
    CtData(t) = D; % Copy most parameters

    CtData(t).parameters.angles = allAngles(t,:);
    CtData(t).parameters.numberImages = Nang;

    decInd = 1 + allAngles(t,:) / angReso; % Index of the desired projection angle
    % Integer (shifted by 1) and fractional part of the random angles
    ind = floor(decInd);
    dec = rem(decInd(:),1); % Forced to column vector!

    % Decimal part gives the opposite amount of weight, i.e. '1.9' is 10% like 1 and 90% like 2
    sinogram(:,:,t) = (1-dec).*OGsino(ind,:) + dec.*OGsino(ind+1,:);
end