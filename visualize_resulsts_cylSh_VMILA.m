% T Heikkilä    2022
% Based on codes by TA Bubba and L Ratti
%% Setup
addpath(genpath('./tools/'));

% Regularization method
p = 1;
transform = 1; % 1:wavelet; 2:cyl.shearlet; 0:identity
%c_alpha = 0.12; % 12
version = 1;

plotsies = true; % Visuals disabled for cluster
redoBregman = false; % Recompute Bregman distances (costly!)

CUDAflag = false;
if CUDAflag
    fprintf('CUDA enabled \n');
else
    fprintf('CUDA disabled \n');
end

% Other
if ~exist('fixed_noise','var')
    fixed_noise = 1;    % 0: delta and alpha reduce as N^{-1}
                        % 1: delta is fixed and alpha reduces as N^{-1/3}
    warning('Experiment parameter missing, set to: %i',fixed_noise)
    fprintf('0 = decreasing, 1 = fixed \n');
end
    
Npxl = size(obj,1);
T = size(obj,3);

% Nsamp = 5; % shorter simulations (max: 30)
Nangsamp = length(allAngles_cell);
% NangMin = 20;
% NangMax = 180;
% 
% if Nangsamp == 1
%     numAngles = NangMin;
% else
%     numAngles = round(linspace(NangMin,NangMax,Nangsamp));
% end


vec = @(x) x(:);
array = @(x,n) reshape(x,n);

% Tile rec_save_cell in a better way
recns = cat(4,rec_save_cell{vec(reshape(1:Nsamp*Nangsamp,Nangsamp,Nsamp)')});

%% Remove some samples
ind = [1:3,6];
BregDist = BregDist(ind,:);
err = err(ind,:);
numAngles = numAngles(ind);
Nangsamp = length(numAngles);

%% Generate transform M
if redoBregman
    fprintf('Recomputing Bregman distances, please wait \n');
    switch transform
        case 1
            xSz = [Npxl,Npxl,T];
            wname = 'db2';
            level = 3;
            wSys = wavedec3(zeros(xSz),level,wname);
            fprintf('Using %s wavelets at level %d \n',wname,level);
        
            getDec  = @(x) x.dec;
            M.times = @(x) getDec(wavedec3(x,level,wname));    % 3D wavelet transform
            M.adj   = @(x) dec2waverec3(x,wSys);    % transpose of 3D wavelet transform
        case 2
            % Listed coarse -> fine (as in theory),
            % algorithm works fine -> coarse (like wavelets)
            decomp=[3 3 4];
            dsize=[16, 16 32];
            level=3 ; % choose level of decomposition ,
            xSz = [Npxl,Npxl,T];
            [shear_f]=setup_cylindrical_shearV2(decomp,dsize,level);

            % Define forward and adjoint for cylindrical shearlets
            % Note the input MUST be reshaped correctly!
    %         M.times = @(x) cylindrical_shearV2(DoPyrDec(x,level),shear_f,level);
            M.times = @(x) cylindrical_shearV2(PyrNDDec_mm(x, 'S', level, 2, @rcos), shear_f,level);
            M.adj   = @(x) cylindrical_shear_adjV2(x,shear_f,level);
        case 0
            M.times = @(x) x;
            M.adj = @(x) x;
    end

    switch p
        case 1
            BregDist_zero = zeros(Nangsamp,Nsamp);
            BregDist = zeros(Nangsamp,Nsamp);
        case 2
            % nothing to do here
        otherwise
            R = @(x) SHpnorm(M.times(x),p); % regularization term             
            signpow = @(x) sign(x).*abs(x).^(p-1);
            gradR = @(x) M.adj(cellfun(signpow,M.times(x),'UniformOutput',false)); % gradient of the regularization term
            bregR = @(u,v) dot(vec(gradR(u)-gradR(v)),vec(u-v));     % Bregman distance
            BregDist = zeros(Nangsamp,Nsamp);
    end
    gradRobj = vec(gradR(obj));
    objDotGradRobj = dot(gradRobj,obj(:));
    for i = 1:Nangsamp
        for j = 1:Nsamp
            rec = rec_save_cell{i,j};
            BregDist(i,j) = bregR(obj,rec) / objDotGradRobj;
            fprintf('.');
        end
        fprintf('\n');
    end
end

%% Comparison plot - Bregman
if plotsies % All plots should be disabled

    % Expected decay
    switch fixed_noise
        case 1
            exp_power = -1/3;
        case 0
            exp_power = -1;
    end
    start = 1;

    indicator = sum(BregDist,2)/Nsamp;
    indicator = indicator';

    % Fitting curve
    yy = log(indicator(start:end)');
    AA = [ones(Nangsamp-start+1,1) log(numAngles(start:end)')];
    cc = AA\yy;
    comparison = exp(AA*cc);
    cc2 = ones(Nangsamp-start+1,1)\(yy-exp_power*log(numAngles(start:end)'));
    expected = exp(AA*[cc2;exp_power]);

    figure(2)
    subplot(1,2,1)
    loglog(numAngles,indicator,'b*-');
    hold on;
    loglog(numAngles(start:end),comparison,'k--');
    loglog(numAngles(start:end),expected,'r-.');

    % Legend
    xlabel('N','interpreter','latex');
    N2 = ['$N^{',num2str(cc(2)),'}$ (fit)'];
    N3 = ['$N^{',num2str(exp_power),'}$ (theo)'];
    leg = legend({'$\mathbf{E} D_R(f_N,f^*)$',N2,N3},'Location','SW');
    set(leg,'Interpreter','latex');
    title(['Breg Dist. Min: ',num2str(min(indicator))]);

    % Comparison plot - error
    indicator = sum(err,2)/Nsamp;
    indicator = indicator';

    % Fitting curve
    yy = log(indicator(start:end)');
    AA = [ones(Nangsamp-start+1,1) log(numAngles(start:end)')];
    cc = AA\yy;
    comparison = exp(AA*cc);
    cc2 = ones(Nangsamp-start+1,1)\(yy-exp_power*log(numAngles(start:end)'));
    expected = exp(AA*[cc2;exp_power]);
    
    set(gca,'fontsize',13);

    figure(2)
    subplot(1,2,2)
    loglog(numAngles,indicator,'b*-');
    hold on;
    loglog(numAngles(start:end),comparison,'k--');
    loglog(numAngles(start:end),expected,'r-.');

    % Legend and axes
    xlabel('N','interpreter','latex');
    N2 = ['$N^{',num2str(cc(2)),'}$ (fit)'];
    N3 = ['$N^{',num2str(exp_power),'}$ (theo)'];
    leg = legend({'$\mathbf{E} D_R(f_N,f^*)$',N2,N3},'Location','SW');
    set(leg,'Interpreter','latex');

    title(['Error. Min: ',num2str(min(indicator))]);

    set(gca,'fontsize',13);
    
end