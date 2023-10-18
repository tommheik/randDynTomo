function [ A ] = create_blkdiag_ct_operator_2d_parallel_astra_cpu( CtData, xDim, yDim, angles )
%CREATE_BLKDIAG_CT_OPERATOR_2D_PARALLEL_ASTRA_CPU Create 2D CT forward model
%    in block diagonal format.
%   A = create_blkdiag_ct_operator_2d_parallel_astra( ctData, xDim, yDim, angles ) computes the 
%   forward model, i.e. X-ray projection operator, for the 2D parallel-beam CT
%   project given in input parameter ''CtData''. The x- and y-dimensions of
%   the CT volume are given by parameters ''xDim'' and ''yDim'', 
%   respectively. The imaging geometry is created using the metadata in 
%   CtData. It is assumed that a flat detector has been used for the X-ray 
%   projection measurements.
%   
%   Each block in the block diagonal operator is created using the 
%   corresponding row of angular values from matrix ''angles''.
%
%   The forward model is an operator that behaves like a matrix, for
%   example in operations like A*x and and A.'*x, but no explicit matrix is
%   actually created.
%
%   Use of this function requires that the ASTRA Tomography Toolbox 
%   (https://www.astra-toolbox.com/) and the Spot Linear-Operator Toolbox 
%   (https://www.cs.ubc.ca/labs/scl/spot/) have been added to the MATLAB 
%   path.
%
%   This function is adapted from the HelTomo Toolbox, which was created 
%   primarily for use with CT data measured in the Industrial Mathematics 
%   Computed Tomography Laboratory at the University of Helsinki.
%
%   T. Heikkilä
%   Created:            9.12.2022
%   Last edited:        18.10.2023
%   
%   Based on codes by 
%   Alexander Meaney, University of Helsinki

% Validate input parameters

if ~isscalar(xDim) || xDim < 1 || floor(xDim) ~= xDim
    error('Parameter ''xDim'' must be a positive integer.');
end

if ~isscalar(yDim) || yDim < 1 || floor(yDim) ~= yDim
    error('Parameter ''yDim'' must be a positive integer.');
end

% Create shorthands for needed variables
if isfield(CtData.parameters,'numDetectors') % v1 field name
    numDetectors    = CtData.parameters.numDetectors;
elseif isfield(CtData.parameters,'numDetectorsPost') % v2 field name
    numDetectors    = CtData.parameters.numDetectorsPost;
end
pixelSize       = CtData.parameters.pixelSize;

% ASTRA uses angles in radians
anglesRad       = deg2rad(angles);

% ASTRA code begins here
fprintf('Creating geometries and data objects in ASTRA... ');

% Create volume geometry, i.e. reconstruction geometry
volumeGeometry = astra_create_vol_geom(xDim, yDim);

numOps = size(angles,1);
if isvector(angles)
    numOps = 1;
    angles = angles(:)'; % Force to row vector
end
Alist = cell(1,numOps);

for i = 1:numOps
    % Create projection geometry
    projectionGeometry = astra_create_proj_geom('parallel', pixelSize, numDetectors, ...
                                                anglesRad(i,:));
    % Create the Spot operator for ASTRA using the GPU.
    Alist{i} = opTomo('strip', projectionGeometry, volumeGeometry);
end

A = blkdiag(Alist{:});

fprintf('done.\n');

% Memory cleanup
astra_mex_data2d('delete', volumeGeometry);
astra_mex_data2d('delete', projectionGeometry);
clearvars -except A

end

