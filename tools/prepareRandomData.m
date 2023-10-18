function [sinogram, CtData, obj] = prepareRandomData(objBig, allAngles)
%PREPARERANDOMDATA Simulates dynamic data given the upsamled phantom
%'objBig' and set of random angles 'allAngles'
% NOTE: this uses parallel beam geometry!
%
%   INPUT
% objBig    Phantom given at twice the spatial resolution N and desired
%           number of time steps T
% allAngles T sequences of random projection angles used for the simulation
%
%   OUTPUT
% sinogram  Sequence of sinograms WITHOUT any noise
% CtData    HelTomo structure for creating the (block diagonal) forward operator
% obj       Downsampled version of objBig
%
% T H   2022



szCell = num2cell(size(objBig)./[2,2,1]);
[N, Ny, T] = szCell{:};
if N ~= Ny
    error('Object should be spatially square! Now it is %d x %d!',N,Ny);
end
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

% CtData structure parameters
binning         = round(1024/N);
pixelSize = 0.1*binning;
% Set desired detector geometry (based on N)
rows            = [36 72 144 288 576 1152];
rows            = min(rows(rows > N));

%% Store parameters in a data structure
% Alexander's code reads these from a text file
CtData = struct;
CtData.type = '2D';
CtData.parameters.angles = [];
CtData.parameters.numberImages = Nang;
CtData.parameters.addedNoise = 0;
CtData.parameters.scanner = 'simulated data';
CtData.parameters.detectorType = 'EID';
CtData.parameters.projectionRows = rows;
CtData.parameters.numDetectors = rows;
CtData.parameters.numDetectorsPost = rows;
CtData.parameters.pixelSizePost = binning*0.1;
CtData.parameters.pixelSize = binning*0.1;
CtData.parameters.binningPost = binning;
CtData.parameters.dataFormat = 'hybrid';


%% Simulate data

% Create volume geometry, i.e. reconstruction geometry
volumeGeometry = astra_create_vol_geom(2*N, 2*N);

sinogram = zeros(Nang,rows,T); % Initialize sinogram

for t = 1:T
    % ASTRA uses radians
    anglesRad = deg2rad(allAngles(t,:));

    % Create projection geometry
    % PixelSize stays the same, there is simply twice as many pixels
    projectionGeometry = astra_create_proj_geom('parallel', pixelSize,2*rows, ...
                                            anglesRad);
    % Create the Spot operator for ASTRA using the CPU
    A = opTomo('strip', projectionGeometry, volumeGeometry);
    
    obj_t = objBig(:,:,t);
    m = reshape(A*obj_t(:),Nang,2*rows); % This sinogram is twice as wide
    sinogram(:,:,t) = 0.25*(m(:,1:2:end-1) + m(:,2:2:end));
    % ASTRA is funky and hence we scale by 0.25 instead of 0.5!
end

% Clear ASTRA memory
astra_mex_data2d('delete', volumeGeometry);
astra_mex_data2d('delete', projectionGeometry);

%% Downsample the true objective
obj = 0.5*(objBig(1:2:end-1,:,:) + objBig(2:2:end,:,:)); % 1st direction
obj = 0.5*(obj(:,1:2:end-1,:) + obj(:,2:2:end,:));       % 2nd direction

