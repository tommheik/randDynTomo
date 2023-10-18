function [x, iter, out] = VMILAn_CTrandomAngles_l1_NNconstraint(gn, R, Phi, beta, opts)
%  VMILA algorithm for shearlets regularized reconstruction problem from CT projections with random angles.
%  version 1.0, T. Bubba, 09/04/2021
%  version 2.0, T. Heikkilä,   28/03/2023 (cylindrical shearlets)
%                  
%   This function solves a tomographic image reconstruction problem by applying 
%   VMILA to the UNCONSTRAINED minimization of the following objection 
%   funcion:
%
%       min 0.5 || R*f - y ||^2 + lam * ||SH(f)||_1 
%        f
%
%   where R is the CT system matrix, y are the observed data, SH 
%   is the wavelet/shearlet matrix, lambda is the regularization parameter 
%   for the regularization term (set in the optional parameters).
%
%
% SYNOPSIS
%   [x, iter, out] = VMILAn_CTrandomAngles_l1(gn, R, Phi[, opts])
%
% MANDATORY INPUT
%   R        (structure of func handle) 
%                           - CT forward operator and its transpose used to  
%                             simulate the syntetic data (sinogram).
%                             If W is a function handle, also the transpose of W is required.
%   gn (y)   (double array) - measured ROI sinogram
%   Phi (SH) (structure of func handle)
%                           - shearlet (or wavelet) operator and its transpose for the 
%                             regularization term                  
%
% OPTIONAL INPUT
%   The following options must be provided as fields for the structure opts.
%   'obj'              - Exact solution, for error calculation (double array)
%   'f0'               - User-provided starting point (double array)
%                        If there is no 'x0' field, a warning message 
%                        occurs.
%   'maxit'            - Maximum number of iterations (integer)
%                        DEFAULT = 1000
%   'beta'             - Beta is the (inverse) of the regularization parameter 
%                        'lam' for the regularization term, placed in front of
%                        the discrepancy term:
%                           min beta/2 * || R*f-y ||^2 + ||\Phi((I-MSK)*Wf + y0))||_1
%                            f
%                        DEFAULT = 1
%   'stopcrit'         - Choice for stopping rule (integer)
%                        1 -> iter > MAXIT
%                        2 -> ||x_k - x_(k-1)|| <= tol*||x_k|| OR iter > MAXIT
%                        3 -> |ObjFun_k - ObjFun_(k-1)| <= tol*|ObjFun_k| OR iter > MAXIT
%                        4 -> ObjFun_k <= tol OR iter > MAXIT
%                        DEFAULT = 1;
%   'tolstop'          - Tolerance used in the stopping criterion
%                        DEFAULT = 1e-4 
%   'eta'              - Eta is the tolerance value for the stopping
%                        criterion of the inner solver. 
%                        If there is no 'eta' field, a warning message 
%                        occurs.
%   'save'             - Directory name: for each iteration saves the 
%                        reconstructed image x_k                       
%   'savefreq'         - Saving frequency (integer or integer array): if it 
%                        is a scalar, then the reconstructed image x_k 
%                        is saved in the directory named in 'save' option.
%                        If it is a positive integer vector of iteration 
%                        counts, then the image and the residual are saved
%                        at those iterations only.                        
%                        If the value/s is/are not integer, then the integer
%                        part is taken. If they are not positive or the 
%                        argument is present but empty, then a warning  
%                        appears and the saving is disabled.
%                        DEFAULT = 1 (save at each iteration)
%
% OUTPUT
%   x (f)              - Reconstructed data
%   iter               - Number of iterations
%
% OPTIONAL OUTPUT
%   The following options must be provided as fields for the structure out.
%   err                - Error value at each iteration.
%                        If 'obj' was not given, then err is an empty matrix.
%   psnr               - PSNR value at each iteration.
%                        If 'obj' was not given, then psnr is an empty matrix.
%   times              - CPU time after each iteration 
%   fv_count           - Computational cost (number of function evaluations).
%   funct              - Value of the objective function at each iteration.
%   InterRec           - Matrix containing successive approximations of x (f) on each column. 
%                        The number of columns depends on savefreq.
%

% Input parameters
if isfield(opts,'x0'), x = opts.x0; else error('Provide initial guess'); end
if isfield(opts,'obj'), obj = opts.obj; flag_obj = 1; else flag_obj = 0; end
% if isfield(opts,'beta'), beta = opts.beta; else beta = 1; end
if isfield(opts,'UpBound'), UpBound = opts.UpBound; boxflag = 1; else boxflag = 0; UpBound = NaN; end
if isfield(opts,'scaling'), scalflag = 1; else scalflag = 0; end
if isfield(opts,'maxit'), maxit = opts.maxit; else maxit = 1000; end
if isfield(opts,'stopcrit'), stopcrit = opts.stopcrit; else stopcrit = 1; end
if isfield(opts,'eta'), eta = opts.eta; else eta = 1e-5; end
if isfield(opts,'tolstop'), tolstop = opts.tolstop; else tolstop = 1e-4; end
if isfield(opts,'save'), save = opts.save; mysave = 1; iter_save  = 0; else mysave = 0; end
if isfield(opts,'savefreq'), savefreq = opts.savefreq; else savefreq = 1e-4; end

if ( ~all(x) )
    initflag = 0;
else
    initflag = 1;
end


if (~isempty(savefreq))
    savefreq = fix( savefreq(:) );
    if ( any(savefreq <= 0) )
        warning('Non-positive saving frequency: saving disabled.')
        mysave = false;
    end
else
    warning('Empty saving frequency: saving disabled.');
    mysave = false;
end


% Steplength default parameters
alpha_min = 1e-5;
alpha_max = 1e5;
gamma     = 1.e-4;
rho       = 0.4;
Malpha    = 3;                          % alphaBB2 memory size
tau       = 0.5;                        % alphaBB alternating parameter
alpha     = 1.3;                        % initial alpha
Valpha    = alpha_max * ones(Malpha,1);
L2        = 1.0e10;                     % bound for scaling matrix


% Vectorization
x_size  = size(x);   % initial guess
x       = x(:);
gn      = gn(:); 

N2 = numel(x);

% output dir handling
if mysave
   [success, mess] = mkdir(save);
   if not(success)
       error('%s: %s',save,mess);
   end
   if not(isempty(mess))
       fprintf('%s\n\n',mess);
   end
%    fh = figure('visible','off');
   InterRec = zeros(prod(x_size),ceil(maxit/savefreq)+1);
end

% Errors and Vector allocation
if flag_obj
    err = zeros(maxit+1,1);
    obj = obj(:);
    obj_norm = sum(obj.*obj);
    e = x - obj;
    err(1) = norm(e(:)) / obj_norm;
end

NrInnerIt = zeros(maxit+1,1);
TimeCost = zeros(maxit+1,1);

% Values of the objective function and gradient at the beginning 
x_tf = R*x;
x_sh = Phi.times(reshape(x,x_size));
d_sh = cellfun(@size,x_sh,'UniformOutput',false);
v1 = cellfun(@zeros,d_sh,'UniformOutput',false);
v2 = zeros(N2,1);
fv = 0.5 * norm(x_tf-gn)^2 + beta * cellnorm(x_sh, 1); % objective function (complete: smooth + nonsmooth)

g = R'*(x_tf-gn);  % gradient of the smooth part

funct = zeros(maxit+1,1);
funct(1) = fv;

% Scaling matrix: adaptive bound
if scalflag
    X_upp_bound = sqrt(1+L2);
    X_low_bound = 1/X_upp_bound;
    
    if initflag == 0
        X = ones(size(x));
    else
        X = x./(R'*x_tf);

        % bounds
        X( X < X_low_bound | isnan(X) ) = X_low_bound;
        X( X > X_upp_bound ) = X_upp_bound;
    end
    % D = 1./X;
else
    X = ones(size(x));
end

iter = 1;
loop = true;
tic
while loop 
    
    fprintf('\n Iteration %d', iter);
    Valpha(1:Malpha-1) = Valpha(2:Malpha);
    
    % Step along the scaled gradient
    z = x - alpha*X.*g;  
        
    % Inexact approximation for p(x^k,h): inner solver (SGP) - composition with a linear operator 
    if iter==1 % Inner solver initialization
        for l = 1:length(d_sh)
            v1{l} = zeros(d_sh{l});
        end
        v2 = zeros(N2,1);
    end
    [y, Delta, v1, v2, cont] = SGP_VMILAn_CTrandomAngles_l1_NNconstraint(X, g, reshape(x, x_size), z, v1, v2, Phi, alpha, eta, beta, boxflag, UpBound);
    NrInnerIt(iter) = cont;  
    fprintf(' Inner Iteration %d', cont);

    d = y - x;
    lambda = 1;
    
    % Backtracking
    fcontinue = 1;

    % Exploiting linearity
    d_tf = R*d;  
    fr = fv;
    
    while fcontinue
        xplus = x + lambda*d;  % new point (try)

        x_tf_try = x_tf + lambda * d_tf;    
        
        fv = 0.5 * norm(x_tf_try-gn)^2 + beta * cellnorm(Phi.times(reshape(xplus,x_size)), 1);
        
        % Sufficient decrease condition
        if ( fv <= fr + gamma * lambda * Delta || lambda<1e-12)
            x = xplus;
            x_tf = x_tf_try;
            sk = lambda*d;
            gtemp = R'*(x_tf-gn);
            yk = gtemp - g;
            g = gtemp;
            fcontinue = 0;
        else
            lambda = lambda * rho;
        end
    end
    
    if (fv >= fr)
        disp('Warning: fv >= fr');
    end
    
    % Update the scaling matrix
    if scalflag
        X = x./(R'*x_tf);
        X_upp_bound = sqrt(1+L2/(iter^2));
        X_low_bound = 1/X_upp_bound;
        X( X < X_low_bound ) = X_low_bound;
        X( X > X_upp_bound ) = X_upp_bound;
        D = 1./X;
    else
        X = ones(size(x));
        D = X;
    end

    sk2 = sk.*D; yk2 = yk.*X;
    bk = dot(sk2(:),yk(:));  ck = dot(yk2(:),sk(:));
    if (bk <= 0)
        alpha1 = min(10*alpha,alpha_max);
    else
        alpha1BB = dot(sk2(:),sk2(:))/bk;
        alpha1 = min(alpha_max, max(alpha_min, alpha1BB));
    end
    if (ck <= 0)
        alpha2 = min(10*alpha,alpha_max);
    else
        alpha2BB = ck/ dot(yk2(:),yk2(:));
        alpha2 = min(alpha_max, max(alpha_min, alpha2BB));
    end
       
    Valpha(Malpha) = alpha2;

    if (iter <= 20)
        alpha = min(Valpha);
    elseif (alpha2/alpha1 < tau)
        alpha = min(Valpha);
        tau = tau*0.9;
    else
        alpha = alpha1;
        tau = tau*1.1;
    end
   
    alpha = double(single(alpha));
    TimeCost(iter+1) = toc;
    
    % Error
    if flag_obj
        e = x - obj;
        err(iter+1) = norm(e(:)) / obj_norm; 
    end
    funct(iter+1) = fv;
    
    % Stop criteria
    switch stopcrit
        case 1
            % Maximum number of iterations
        case 2
            % || x_k - x_k-1 ||_2 / || x_k ||_2
            crit_value = norm(sk)/norm(x);
            loop = crit_value > tolstop;
        case 3
            % | f_k - f_k-1 | / |f_k|
            crit_value = abs(fv - funct(iter))/abs(fv); 
            loop = crit_value > tolstop;
        case 4
            % f_k < tol
            crit_value = fv * scaling^2;
            discr = 0.5*x_size(1,1)*x_size(1,2);
            loop = crit_value > discr;
     end
    
     if iter >= maxit
         loop = false;
     end
     if mysave
         if ( length(savefreq) > 1 && any(iter-1 == savefreq) ) || ( length(savefreq) == 1 && ~rem(iter-1, savefreq) )
             iter_save = iter_save + 1;
             InterRec(:,iter_save) = x;
             savedLastIter = 1;
         else
             savedLastIter = 0;
         end
     end
       
    iter = iter + 1;        
end

% save final images
if mysave
    if ( ~savedLastIter ) 
        InterRec(:,iter_save+1) = x;
    end
    out.InterRec = InterRec; 
end


%% Output
x = reshape(x, x_size);

iter = iter - 1;

if flag_obj
    out.err = err(1:iter);
end
out.funct = funct(1:iter);
out.NrInnerIt = NrInnerIt(1:iter-1);
out.TimeCost = TimeCost(1:iter);

end






