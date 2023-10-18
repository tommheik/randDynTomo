% version 1.0, T. Bubba, 10 August 2021
% version 2.0, T H,   18 November 2022 (Updated 18.6.2023)

close all
clear 

dataset = 'Stempo'; % 'Gel'; % 'Cartoon'

dateString = string(datetime('now','TimeZone','local','Format','yyMMdd-HHmmss'));
fprintf('Begin at: %s \n',dateString);

fprintf('Reconstructing dataset: %s \n', dataset);

% Only the old ASTRA toolbox works
addpath(genpath('../PDFP4D/tools/astra-toolbox'))
addpath(genpath('./tools/'));

% Regularization method
p = 1.5; % 1.5;
transform = 1; % 1:wavelet; 2:cyl.shearlet; 0:identity
version = 1;

plotsies = false; % Visuals disabled for cluster

CUDAflag = false;
if CUDAflag
    fprintf('CUDA enabled \n');
else
    fprintf('CUDA disabled \n');
end

switch dataset
    case 'Cartoon'
        if transform == 2 % cylindrical shearlet
            c_alpha = 0.03;
        else % Wavelet
            c_alpha = 0.12;
        end

        obj = load('./phantoms/cartoonPhantom_256x256x32.mat');
        objBig = obj.obj; % Upsampled phantom
        obj = 0.5*(objBig(1:2:end-1,:,:) + objBig(2:2:end,:,:));
        obj = 0.5*(obj(:,1:2:end-1,:) + obj(:,2:2:end,:));

    case 'Gel'
        if transform == 2 % cylindrical shearlet
            c_alpha = 0.005;
        else % Wavelet
            c_alpha = 0.02;
        end

        T = 16;
        % Matlab loads the variables as structs containing the variables
        CtDataOriginal = load('./gelPhantom/GelPhantomData_b4.mat','GelPhantomData_b4');

        % Choose which reconstruction is used as 'ground truth' for computing the Bregman distances
        objFile = 'CIL_tv_256x256x17_scaled.mat'; % 'fbp_256x256x17.mat'
        fprintf('Using file: %s as ground truth \n',objFile)
        obj = load(['./gelPhantom/', objFile],'obj');

        CtDataOriginal = CtDataOriginal.GelPhantomData_b4(1:T);
        obj = rot90(obj.obj(:,:,1:T));

    case 'Stempo'
        if transform == 2 % cylindrical shearlet
            c_alpha = 0.005;
        else % Wavelet
            c_alpha = 0.1;
        end

        T = 16;
        % Matlab loads the variables as structs containing the variables
        CtDataOriginal = load('./stempo/stempo_seq8x45_2d_b8.mat');
        CtDataOriginal = CtDataOriginal.CtData;

        % Stupid things need to be done
        CtDataOriginal.parameters.numDetectors = CtDataOriginal.parameters.numDetectorsPost;
        CtDataOriginal.parameters.pixelSize = CtDataOriginal.parameters.pixelSizePost;
        CtDataOriginal.parameters.effectivePixelSize = CtDataOriginal.parameters.effectivePixelSizePost;

        % Choose which reconstruction is used as 'ground truth' for computing the Bregman distances
        objFile = 'stempo_ground_truth_2d_b4.mat'; % 'fbp_256x256x17.mat'
        fprintf('Using file: %s as ground truth \n',objFile)
        obj = load(['./stempo/', objFile],'obj');
        % We only pick 16 time steps for reference
        objTimeSteps = [2 25 49 73 97 120 144 169 193 217 240 264 287 311 333 328];
        obj = obj.obj(:,:,objTimeSteps);

        fprintf('Resizing obj to 280x280\n');
        obj = imresize(obj,[280,280]);

end

% Other
fixed_noise = 0;    % 0: delta and alpha reduce as N^{-1}
                    % 1: delta is fixed and alpha reduces as N^{-1/3}
Npxl = size(obj,1);
T = size(obj,3);

Nsamp = 5; % shorter simulations (max: 30)

%% Choose number of angles N
switch dataset
    case 'Stempo'
        % Not so many angles available for Stempo data
        Nangsamp = 6;
        NangMin = 8;
        NangMax = 50;
    otherwise
        % Shearlet: 16:240 with 8
        % Wavelet: 16:520 with 10 (not that good actually)
        Nangsamp = 8;
        NangMin = 16;
        NangMax = 240;
end

if Nangsamp == 1
    numAngles = NangMin;
else
    % Even spacing in log-log plots
    numAngles = round(exp(linspace(log(NangMin),log(NangMax),Nangsamp)));
    fprintf('Number of angles used:\n')
    disp(numAngles)
end


vec = @(x) x(:);
array = @(x,n) reshape(x,n);


%% Generate transform M
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
        decomp=[3 3 4]; % [2 3 3 4];
        dsize=[16, 16, 32]; % [8 16, 16 32];
        level=3; % 4; % choose level of decomposition ,
        xSz = [Npxl,Npxl,T];
        [shear_f]=setup_cylindrical_shearV2(decomp,dsize,level);
        fprintf('Using cylindrical shearlets at level %d \n', level);
        
        % Define forward and adjoint for cylindrical shearlets
        % Note the input MUST be reshaped correctly!
%         M.times = @(x) cylindrical_shearV2(DoPyrDec(x,level),shear_f,level);
        M.times = @(x) cylindrical_shearV2(PyrNDDec_mm(x, 'S', level, 2, @rcos), shear_f,level);
        M.adj   = @(x) cylindrical_shear_adjV2(x,shear_f,level);
    case 0
        M.times = @(x) x;
        M.adj = @(x) x;
end

%% VMILA optional parameters
opts = struct();
opts.x0       = zeros(xSz);  
normobj = norm(obj(:));          % norm of the phantom
opts.obj      = obj;
opts.tolstop  = 1e-5;
opts.stopcrit = 2;
opts.scaling  = true;
opts.savefreq = 5;
% opts.maxit = 3;

%% Allocate vectors to store values for final plots
err = zeros(Nangsamp,Nsamp);
%fprintf('p = %s \n', string(sym(p)))
switch p
    case 1
        Mobj = M.times(obj); % Transform of the ground truth f
        fprintf('Precomputing ||M obj||')
        MobjNorm = cellnorm(Mobj,1);

        bregR = @(Mu,Mv) bregL1(Mu, Mv);
        
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
DataDiscr = zeros(Nangsamp,Nsamp);
rec_save_cell = cell(Nangsamp,Nsamp); %zeros(Npxl^2,Nangsamp);

allAngles_cell = cell(Nangsamp,1);


%% Angle dependent parameters

switch fixed_noise % Note mMax is missing from c_delta!
    case 1
        % Fixed noise level delta, alpha decreases as angles increase
        c_delta = 0.03; % 0.03 good, 0.01 was not good
    case 0
        % Both alpha and delta decrease as angles increase
        c_delta = 0.035*numAngles(1);
end

switch fixed_noise
    case 1
        alpha = c_alpha*numAngles.^(-1/3);
        delta = c_delta*ones(1,Nangsamp);
    case 0
        alpha = c_alpha*numAngles.^(-1);
        delta = c_delta.*numAngles.^(-1);
end

for i = 1:Nangsamp
    Nang = numAngles(i); % Number of projection angles per time step
    fprintf('Number of angles %d (%d out of %d): \n',Nang,i,Nangsamp)

    %% Generate data, operator and apply noise
    rngSeed = 42; % Used to make everything 'random' repeatable (and iterable with different seed)

    fprintf('Fixing (angular) randomness with seed: %d! \n',rngSeed)
    rng(rngSeed,'twister')

    % Generate random angles for each repetition and each time step
    switch dataset
        case 'Stempo'
            % Stempo uses different setup, projections only go 180 degrees
            Tshift = [180*(0:1:T-2),(T-1)*180 - 8]; % Last time step is shifted only 172 degrees!
            allAngles = 180*rand(Nsamp, T, Nang) + Tshift;
        otherwise
            % other data is simpler
            allAngles = 360*rand(Nsamp, T, Nang);
    end
    allAngles_cell{i} = allAngles; % Store angles
    
    for k = 1:Nsamp % Repeat same number of angles but different actual angles and realizations of noise
        fprintf('Sample %d (%d angles): \t',k,Nang)
        
        Angles = squeeze(allAngles(k,:,:)); % Unique angles
        fprintf('Angles used:\n      %1.2f, %1.2f, ..., %1.2f (deg)\n      ...\n      %1.2f, %1.2f, ..., %1.2f \n',...
            Angles(1,1), Angles(1,2), Angles(1,end), Angles(end,1), Angles(end,2), Angles(end,end));
        
        switch dataset
            case 'Cartoon'
                % Generate noise-free data for these projection angles
                [m, CtData] = prepareRandomData(objBig,Angles);
                CtData.parameters.angles = Angles(1,:); % For generating one operator for t=1

                % Compute norm of A now for efficiency
                Asmall = create_ct_operator_2d_parallel_astra(CtData,Npxl,Npxl);
                normA = 1/normest(Asmall); % Spectral norm (p=2) of a block diagonal matrix = the spectral norm of one of the blocks
                % Important: normest uses a lot of random vectors!

                % Block diagonal operator with the random angles
                if CUDAflag
                    error('Parallel beam CUDA implementation missing!')
                    A = create_blkdiag_ct_operator_2d_fan_astra_cuda(CtData,Npxl,Npxl,Angles); 
                else
                    A = normA*create_blkdiag_ct_operator_2d_parallel_astra_cpu(CtData,Npxl,Npxl,Angles); 
                    % BlockDiag operator weighted by 1/norm(A) to normalize
                    % A = opBlockDiag(normA*ones(1,T),Asmall);
                end

            case 'Gel'
                [m, CtData] = extrapolateCtData(CtDataOriginal, Angles);

                % Compute norm of A now for efficiency
                Asmall = create_ct_operator_2d_fan_astra(CtData(1),Npxl,Npxl);
                normA = 1/normest(Asmall); % Spectral norm (p=2) of a block diagonal matrix = the spectral norm of one of the blocks
                % Important: normest uses a lot of random vectors!

                % Block diagonal operator with the random angles
                if CUDAflag
                    A = create_blkdiag_ct_operator_2d_fan_astra_cuda(CtData(1),Npxl,Npxl,Angles); 
                else
                    A = normA*create_blkdiag_ct_operator_2d_fan_astra_cpu(CtData(1),Npxl,Npxl,Angles); 
                end

            case 'Stempo'
                [m, CtData] = extrapolateStempoCtData(CtDataOriginal, Angles);

                % Compute norm of A now for efficiency
                Asmall = create_ct_operator_2d_fan_astra(CtData(1),Npxl,Npxl);
                normA = 1/normest(Asmall); % Spectral norm (p=2) of a block diagonal matrix = the spectral norm of one of the blocks
                % Important: normest uses a lot of random vectors!

                % Block diagonal operator with the random angles
                if CUDAflag
                    A = create_blkdiag_ct_operator_2d_fan_astra_cuda(CtData(1),Npxl,Npxl,Angles); 
                else
                    A = normA*create_blkdiag_ct_operator_2d_fan_astra_cpu(CtData(1),Npxl,Npxl,Angles); 
                end
        end

        % Normalize
        m = m * normA;
        mMax = max(m(:));
        
        dateString = string(datetime('now','TimeZone','local','Format','yyMMdd-HHmmss'));
        fprintf('Current time: %s \n',dateString);
        
        fprintf('Reset (noise) randomness with seed: %d! \n',k*rngSeed) % Needs to change for the noise to be different
        rng(k*rngSeed,'twister')
        noise = randn(size(m));
        
        % Noisy measurements (notice how mMax is included here!)
        meas = m + delta(i)*noise*mMax;
        
        %% Begin outer iteration
        
        % diary(savename)
        switch p
            case 1
                [rec,~,out] = VMILAn_CTrandomAngles_l1_NNconstraint(meas, A, M, alpha(i), opts);
                BregDist(i,k) = bregL1(Mobj, M.times(rec)) / MobjNorm;
%                 BregDist_zero(i,k) = SymbregR2(f(:),rec(:),M)/(SymbregR2(f(:),0*rec(:),M));
            case 2
                error('Not implemented yet!')
                J = @(f) 1/2*norm(A.times(f)-meas)^2 + alpha(i)/2*norm(f)^2;
                gradJ = @(f) A.adj(A.times(f)-meas) + alpha(i)*f;
                [rec,~,~,~,~] = SGP_generic(J, gradJ, opts.x0);
                BregDist(i,k) = 1/2*norm(rec(:)-f(:))^2/normobj^2;
            otherwise
                [rec,~,out] = VMILAn_CTrandomAngles_lp_NNconstraint(meas, A, M, p, alpha(i), opts);
                BregDist(i,k) = bregR(obj,rec) / dot(vec(gradR(obj)),obj(:));%(normgradf*normf);
        end
        err(i,k) = norm(rec(:)-obj(:))^2/normobj^2;
        DataDiscr(i,k) = norm(A*(vec(rec))-vec(meas))^2 / norm(vec(meas))^2;
        rec_save_cell{i,k} = rec;
        fprintf('Finished with reconstruction!\n')
        
        if k == Nsamp && plotsies
            rec_save(:,i) = rec(:);
            figure(1)
            clf
            subplot(1,2,1)
            imshow(obj(:,:,round(T/2)))
            subplot(1,2,2)
            truesize([300 300])
            imshow(rec(:,:,round(T/2)),[])
            title(['Angles: ',num2str(numAngles(i)),'(',num2str(i),'/',num2str(Nangsamp),'), Sample: ',num2str(k),'/',num2str(Nsamp)]);
            drawnow
            pause(0.3)
        end
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
    loglog(numAngles(start:end),expected,'k-.');

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

    figure(2)
    subplot(1,2,2)
    loglog(numAngles,indicator,'b*-');
    hold on;
    loglog(numAngles(start:end),comparison,'k--');
    loglog(numAngles(start:end),expected,'k-.');

    % Legend and axes
    xlabel('N','interpreter','latex');
    N2 = ['$N^{',num2str(cc(2)),'}$ (fit)'];
    N3 = ['$N^{',num2str(exp_power),'}$ (theo)'];
    leg = legend({'$\mathbf{E} D_R(f_N,f^*)$',N2,N3},'Location','SW');
    set(leg,'Interpreter','latex');

    title(['Error. Min: ',num2str(min(indicator))]);

    set(gca,'fontsize',12);
end

%% Save results

dateString = string(datetime('now','TimeZone','local','Format','yyMMdd-HHmmss'));
fprintf('Finished at: %s \n\n',dateString);

% Save name
savename = ['Results_',dataset];
switch fixed_noise
    case 0
        savename = [savename,'_decreasing'];
    case 1
        savename = [savename,'_fixed'];
end
switch p
    case 1
        savename = [savename,'_p1'];
    case 1.1
        savename = [savename,'_p11'];
    case 1.01
        savename = [savename,'_p101'];    
    case 3/2
        savename = [savename,'_p32'];
    case 4/3
        savename = [savename,'_p43'];
    case 2
        savename = [savename,'_p2'];
end
savename = [savename,'_img',num2str(Npxl),'_Nsamp',num2str(Nsamp,'%03d')];
switch transform
    case 1
        savename = [savename,'_WaveletConstrained'];
		methdFldr = 'wavelet/';
    case 2
        savename = [savename,'_ShearletConstrained'];
		methdFldr = 'shearlet/';
    case 0
        savename = [savename,'_TikhonovConstrained'];
		methdFldr = 'tikhonov/';
end
savename = strcat(savename,'_v',num2str(version),'_Cluster_',dateString,'.mat');

fprintf('Saving as: %s ',savename)

save(fullfile('./results/',methdFldr,savename),'BregDist','DataDiscr','err','opts','fixed_noise',...
    'numAngles','c_alpha','c_delta','obj','Nsamp','rec_save_cell','allAngles_cell','dateString', 'objFile')

fprintf('   complete!\n\n')

