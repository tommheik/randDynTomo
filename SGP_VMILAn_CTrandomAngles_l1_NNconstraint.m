function [prim, Delta, v1, v2, cont] = SGP_VMILAn_CTrandomAngles_l1_NNconstraint(X, grad_smooth, x, z, v1, v2, Phi, Alpha_k, eta, Beta_reg, boxflag, UpBound)
% SGP_VMILA - SGP algorithm adapted to solve the dual problem from eta-approx 
%   of VMILA algorithm:   
%
%   max  -(2^alpha_k)^(-1)*||alpha_k*R^T*Phi^T*v - z_k||_2^2 - f_1(x_k) - 0.5*alpha_k ||Grad(f_0(x_k))||_2^2 + (2^alpha_k)^(-1)*||z_k||_2^2
%     s.t. ||v||_inf <= Beta_reg 
%
%   where R is the system matrix for CT system, 
%   v is the dual variable (primal: y = z_k - alpha_k*W^T*Phi^T*v), 
%   Phi is the wavelet/shearlet matrix, f_1(x) = norm(Phi(x), 1),
%   z_k = x_k - alpha_k*Grad(f_0(x_k)), alpha_k is the steplength and the    
%   feasible region is the unit ball in the infinity norm.
%   In the objective function above the scaling matrix is missing, but it
%   is present in the implementation in this code.
%
% SYNOPSIS
%   [primal, Delta, dual, cont] = SGP_VMILAn_CTrandomAngles_l1(X, grad_smooth, x, z, dual, Phi, alpha, eta, beta, boxflag, UpBound)
%
% MANDATORY INPUT												  
%   X           (double array)
%                           - (diagonal) scaling matrix for variable metric
%   grad_smooth (double array)
%                           - gradient of the smooth part f_0 of the original
%                             objective function ( f_0(x) + f_1(x) )
%   x           (double array)
%                           -  reconstructed data
%   z           (double array)
%                           -  z_k = x_k - alpha_k*Grad(f_0(x_k))
%   v           (double array)
%                           - inizial guess for dual variable (v)
%   Phi         (double matrix or structure with func handle)
%                           - wavelet or shearlet matrix 
%   alpha       (double array)
%                           -  steplength
%   eta         (double array)
%                           -  parameter for stopping criterion
%
% OUTPUT
%   prim               - primal variable (y in the article)
%                        y = z_k - alpha_k*W^T*(I-MSK)^T*Phi^T*v
%   cont               - Number of iterations
%   Delta              - 
%                        Delta = 
%   v                  - dual variable (starting approx for outer loop)
%                        
%


% SGP default parameters 
%   'M'                - Nonmonotone lineasearch memory (positive integer)
%                        If M = 1 the algorithm is monotone
%                        DEFAULT = 1 (monotone)
%   'GAMMA'            - Linesearch sufficient-decrease parameter (double)
%                        DEFAULT = 1e-4
%   'BETA'             - Linesearch backtracking parameter (double)
%                        DEFAULT = 0.4
%   'ALPHA_MIN'        - Lower bound for Barzilai-Borwein' steplength (double>0)
%                        DEFAULT = 1e-5
%   'ALPHA_MAX'        - upper bound for Barzilai-Borwein' steplength (double>0)
%                        DEFAULT = 1e5
%   'MALPHA'           - Memory length for alphaBB2 (positive integer)
%                        DEFAULT = 3
%   'TAUALPHA'         - Alternating parameter for Barzilai-Borwein' steplength
%                        (double)
%                        DEFAULT = 0.5
%   'INITALPHA'        - Initial value for Barzilai-Borwein' steplength (double)
%                        DEFAULT = 1.3 
gamma = 1e-4;          % for sufficient decrease
beta = 0.4;            % backtracking parameter
M = 1;                 % memory in obj. function value (if M = 1 monotone)
alpha_min = 1e-5;      % alpha lower bound
alpha_max = 1e5;	   % alpha upper bound
Malpha = 3;            % alfaBB1 memory
tau = 0.5;             % alternating parameter
initalpha = 1.3;       % initial alpha
maxinnerit = 20;

xSz = size(x);

% Vectorization
x = x(:);
v2 = v2(:);
z = z(:);

% some computations needed only once
cont   = 1;                              % iteration counter
alpha  = initalpha;                      % initial alpha
Valpha = alpha_max * ones(Malpha,1);     % memory buffer for alpha
Fold   = -1e30 * ones(M, 1);             % memory buffer for obj. func.
D = 1./X;

% projection of the initial point (this is how we enforce g^*)
v1 = SH2feas(v1, Beta_reg);
v2 = min(0, v2);

% Compute the primal variable
prim = z - Alpha_k * X.* (vec(Phi.adj(v1)) + v2);

% objective function (fv) and gradient (g) value (scaled version) 
fun_discr = 1/(2*Alpha_k) * (prim'* (D.*prim));  % This is the only piece of the objective function that it is NOT constant wrt the variable v!
x_sh = Phi.times(reshape(x,xSz));
fun_fix = - 1/(2*Alpha_k) * (z'*(D.*z)) + (0.5*Alpha_k) * (grad_smooth'*(X.*grad_smooth)) + Beta_reg * cellnorm(x_sh, 1);
if isnan(fun_fix)
    error('NaN value spotted!')
end
fv = fun_discr + fun_fix; 

% GRADIENT: A(alpha*Dk*(AT*v1 + v2) - z) 
g1 = Phi.times(reshape(-prim,xSz));
g2 = -prim;
prim( prim < 0 ) = 0;
prim = max(0, prim);
if boxflag 
    prim = min(UpBound, prim);
end

% main loop
loop = true;

while loop
    % store alpha and objective function values
    Valpha(1:Malpha-1) = Valpha(2:Malpha);
    Fold(1:M-1) = Fold(2:M);
    Fold(M) = fv;

    % compute descent direction (d)
    y1 = SHsum(v1,g1,1,-alpha); % y1 = v1 - alpha*g1;
    y2 = v2 - alpha*g2;
    
    % projection onto the feasible set
    y1 = SH2feas(y1, Beta_reg);
    y2( y2 > 0 ) = 0;

    d1 = SHsum(y1, v1, 1, -1); % d1 = y1 - v1;
    d2 = y2 - v2;

    % backtracking loop for linesearch
    gd = dot(d2,g2);
    for l = 1:length(d1)
        gd = gd  + dot(d1{l}(:), g1{l}(:)); % + dot(d1,g1)
    end
    lam = 1;
    
    fcontinue = 1;    
    fr = max(Fold);

    while fcontinue
        vplus1 = SHsum(v1, d1, 1, lam);     % dual variable: try update
        vplus2 = v2 + lam*d2;
        prim_try = z - Alpha_k * X.* (vec(Phi.adj(vplus1)) + vplus2); % primal variable: try update
        
        % OBJECTIVE FUNCTION
        fv = 1/(2*Alpha_k) * (prim_try'*(D.*prim_try)) + fun_fix;
        
        if ( fv <= fr + gamma * lam * gd || lam < 1e-12)            
            v1 = vplus1;
            v2 = vplus2;
            sk2 = lam*d2; % sk = lam*[d1; d2];
            gtemp1 = Phi.times(reshape(-prim_try, xSz));
            gtemp2 = -prim_try;
            
            y2 = gtemp2 - g2; % yk = [gtemp1; gtemp2] - [g1; g2];
            y1 = SHsum(gtemp1, g1, 1, -1);
            g1 = gtemp1;
            g2 = gtemp2;
            fcontinue = 0;
        else
            lam = lam * beta;
        end
    end
    prim = prim_try;
    prim = max(prim, 0);
    if boxflag
        prim = min(prim, UpBound);
    end

    % update the steplength 
    bk = dot(sk2(:),y2(:));  % dot(sk1, y1);
    for l = 1:length(d1)
        bk = bk + dot(lam*d1{l}(:), y1{l}(:));
    end
    if (bk <= 0)
        alpha1 = min(10*alpha,alpha_max);
        alpha2 = min(10*alpha,alpha_max);
    else
        alpha1BB = (norm(sk2(:)).^2 + lam^2.*cellnorm(d1,2).^2)/bk; % dot(sk(:),sk(:))/bk;
        alpha1 = min(alpha_max, max(alpha_min, alpha1BB));
        alpha2BB = bk / (norm(y2(:)).^2 + cellnorm(y1,2).^2); % bk/ dot(yk(:),yk(:));
        alpha2 = min(alpha_max, max(alpha_min, alpha2BB));
    end
       
    Valpha(Malpha) = alpha2;

    if (cont <= 20)
        alpha = min(Valpha);
    elseif (alpha2/alpha1 < tau)
        alpha = min(Valpha);
        tau = tau*0.9;
    else
        alpha = alpha1;
       
        tau = tau*1.1;
    end
    
    alpha = double(single(alpha));

    cont = cont + 1;


    % stopping criterion
    dif = prim - x;
    Delta = dot(dif(:),grad_smooth(:)) + 1/(2*Alpha_k)*(dif'*(D.*dif)) + Beta_reg * cellnorm(Phi.times(reshape(prim,xSz)), 1) - Beta_reg * cellnorm(Phi.times(reshape(x,xSz)), 1);
    Psi = -fv;  % Psi is the objective function of the dual problem (ie, the opposite of the objective funtion fv that arise
                % from the min problem equivalent to max dual problem with Psi): don't mistake the optimization problem with the stopping criterion
    if  Delta <= eta*Psi || cont >= maxinnerit
        loop = false;
    end
    
    fprintf('\n         lambda = %e, Delta = %e, eta*Psi=%e\n', lam, Delta, eta*Psi)
    
    
end

cont = cont - 1;


return