function [prim, Delta, v1, v2, cont] = SGP_VMILAn_CTrandomAngles_lp_NNconstraint(X, grad_smooth, x, z, v1, v2, Phi, Alpha_k, eta, Beta_reg, p, boxflag, UpBound, xSz)
% SGP_VMILA - SGP algorithm adapted to solve the dual problem from eta-approx 
%   of VMILA algorithm:   
%
%   max  -(2^alpha_k)^(-1)*||alpha_k*R^T*Phi^T*v - z_k||_2^2 - g^*(v) - f_1(x_k) - 0.5*alpha_k ||Grad(f_0(x_k))||_2^2 + (2^alpha_k)^(-1)*||z_k||_2^2
%
%   where R is the system matrix for CT system, 
%   v is the dual variable (primal: y = z_k - alpha_k*W^T*Phi^T*v), 
%   Phi is the wavelet/shearlet matrix, f_1(x) = 1/p * norm(Phi(x), p)^p,
%   z_k = x_k - alpha_k*Grad(f_0(x_k)), alpha_k is the steplength and the    
%   feasible region is the unit ball in the q-norm where 1/p+1/q = 1.
%   In the objective function above the scaling matrix is missing, but it
%   is present in the implementation in this code.
%
% SYNOPSIS
%   [primal, Delta, dual, cont] = SGP_VMILAn_CTrandomAngles_l1(R, X, grad_smooth, x, z, dual, Phi, alpha, eta, beta, boxflag, UpBound)
%
% MANDATORY INPUT
%   R           (double array or structure with func handle) 
%                           - system matrix for CT system 
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
%   xSz         (double array)
%                           -  size of array x
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
maxinnerit = 50;


% Compute q
if p == 3/2
    q = 3;
elseif p == 4/3
    q = 4;
elseif p == 1.1
    q = 11;
elseif p == 1.01
    q = 101;
else
    error('Pair (p,q) non implemented.')
end

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

% projection of the initial point
v2( v2 > 0 ) = 0;

% Compute the primal variable
prim = z - Alpha_k * X.* (vec(Phi.adj(v1)) + v2);
% fprintf('max(prim) = %5.4f \n ',max(prim(:)))

% objective function (fv) and gradient (g) value (scaled version) 
fun_discr = 1/(2*Alpha_k) * (prim'* (D.*prim)) + (Beta_reg)^(1-q) * SHpnorm(v1,q);  % This is the only piece of the objective function that it is NOT constant wrt the variable v!
% fun_discr = 1/(2*Alpha_k) * (prim'* (D.*prim)) + (Beta_reg)/q * norm(v1,q)^q;  % This is the only piece of the objective function that it is NOT constant wrt the variable v!
x_sh = Phi.times(reshape(x,xSz));
x_sh_pnorm = SHpnorm(x_sh,p);
fun_fix = - 1/(2*Alpha_k) * (z'*(D.*z)) + (0.5*Alpha_k) * (grad_smooth'*(X.*grad_smooth)) + Beta_reg * x_sh_pnorm;
fv = fun_discr + fun_fix; 
% fprintf('funct obj = %5.4f \n ',fv)

% GRADIENT: A(alpha*Dk*(AT*v1 + v2) - z) 
signpowq = @(x) sign(x).*abs(x).^(q-1);           % p-power operation from paper
SHsum = @(u,v,a,b) a*u + b*v; % Scaled sum on shearlet coefficient domain (i.e. cell arrays)
prim_sh = Phi.times(reshape(-prim,xSz));
g1 =  cellfun(@(u,v) SHsum(u,signpowq(v),1,Beta_reg^(1-q)),prim_sh,v1,'UniformOutput',false);
% g1 = Phi.times(-prim) + Beta_reg *signpowq(v1,q-1);
% fprintf('max(grad) = %5.4f \n ',max(g1(:)))
g2 = -prim;
prim( prim < 0 ) = 0;
if boxflag 
    prim( prim > UpBound ) = UpBound;
end


% projection of the initial point (use Frank-Wolfe method: need gradient)
% gam = 1;
% tmp = signpowFW(g1,qp);
% tmp = tmp./norm(tmp,p)*Beta_reg;    % tmp should be s.t. |tmp|_p = Beta_reg
% v1 = (1 - gam) * v1 + gam * tmp;

% main loop
loop = true;

while loop
    % store alpha and objective function values
    Valpha(1:Malpha-1) = Valpha(2:Malpha);
    Fold(1:M-1) = Fold(2:M);
    Fold(M) = fv;

    % compute descent direction (d)
    % y1 = v1 - alpha*g1;
    y1 = cellfun(@(u,v) SHsum(u,v,1,-alpha), v1, g1, 'UniformOutput', false);
    % fprintf('max(y1) = %5.4f \n ',max(y1(:)))
    y2 = v2 - alpha*g2;
    
    % projection onto the feasible set
    y2( y2 > 0 ) = 0;
    
    % projection onto the feasible set
%     gam = 2/(cont+1);
%     tmp = signpowFW(g1,qp);
%     tmp = tmp./norm(tmp,p)*Beta_reg;    % tmp should be s.t. |tmp|_p = Beta_reg
%     y1 = (1 - gam) * y1 + gam * tmp;

    % d1 = y1 - v1;
    d1 = cellfun(@(u,v) SHsum(u,v,1,-1), y1, v1, 'UniformOutput', false);
    % fprintf('max(d1) = %5.4f \n ',max(d1(:)))
    d2 = y2 - v2;

    % backtracking loop for linesearch
    d1dotg1 = sum(cellfun(@(u,v) dot(vec(u),vec(v)), d1, g1));
    gd = d1dotg1 + dot(d2,g2);
    % fprintf('max(gd) = %5.4f \n ',max(gd))
    % fprintf('max(dot(d2,g2)) = %5.4f \n ',max(dot(d2,g2)))
    lam = 1;
    
    fcontinue = 1;    
    fr = max(Fold);

    while fcontinue
        % vplus1 = v1 + lam*d1;     % dual variable: try update (400/9)
        vplus1 = cellfun(@(u,v) SHsum(u,v,1,lam), v1, d1, 'UniformOutput', false);
        vplus2 = v2 + lam*d2;
        prim_try = z - Alpha_k * X.* (vec(Phi.adj(vplus1)) + vplus2); % primal variable: try update 400/9 * 
        % fprintf('max(prim_try) = %5.4f \n ',max(prim_try(:)))
        
        % OBJECTIVE FUNCTION
        fv1 = 1/(2*Alpha_k) * (prim_try'*(D.*prim_try));
        % fprintf('fv1 = %5.4f \n ', fv1)
%         fv2 = ((Beta_reg)^(1-q)/q * norm(vplus1,q)^q); % ./373248
        fv2 = (Beta_reg)^(1-q) * SHpnorm(vplus1,q);
        % fv2 = Beta_reg/q * norm(vplus1,q)^q;
        % fprintf('fv2 = %5.4f \n ', fv2)
        fv3 = fun_fix;
        % fprintf('fv3 = %5.4f \n ', fv3)
        fv = fv1 + fv2 + fv3;
        % fprintf('fv = %5.4f \n ', fv)
        
        if ( fv <= fr + gamma * lam * gd || lam < 1e-12)            
            v1 = vplus1; clear vplus1;
            v2 = vplus2; clear vplus2;
            % sk = lam*[d1; d2];
            sk2 = lam*d2;
            prim_try_sh = Phi.times(reshape(-prim_try,xSz));
            gtemp1 = cellfun(@(u,v) SHsum(u,signpowq(v),1,Beta_reg^(1-q)),prim_try_sh,v1,'UniformOutput',false);
            % gtemp1 = Phi.times(-prim_try) + Beta_reg *signpowq(v1,q-1);
            % fprintf('max(gtemp1) = %5.4f \n ',max(gtemp1(:)))
%             gtemp2 = -prim_try;
            
            % yk = [gtemp1; gtemp2] - [g1; g2];
            yk1 = cellfun(@(u,v) SHsum(u,v,1,-1), gtemp1, g1, 'UniformOutput', false);
            yk2 = -prim_try - g2;
            g1 = gtemp1; clear gtemp1;
            g2 = -prim_try;
            fcontinue = 0;
        else
            lam = lam * beta;
        end
    end
    prim = prim_try;
    prim( prim < 0 ) = 0;
    if boxflag
        prim( prim > UpBound ) = UpBound;
    end

    % update the steplength 
    % bk = dot(sk, yk);  
    bk = dot(sk2,yk2) + sum(cellfun(@(u,v) dot(lam*vec(u),vec(v)), d1, yk1));
    if (bk <= 0)
        alpha1 = min(10*alpha,alpha_max);
        alpha2 = min(10*alpha,alpha_max);
    else
        skNorm = norm(sk2)^2 + sum(cellfun(@(x) norm(lam*x(:))^2, d1));
        ykNorm = norm(yk2)^2 + sum(cellfun(@(x) norm(x(:))^2, yk1));
        alpha1BB = skNorm/bk;
        alpha1 = min(alpha_max, max(alpha_min, alpha1BB));
        alpha2BB = bk/ ykNorm;
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
%     Delta = dot(dif(:),grad_smooth(:)) + 1/(2*Alpha_k)*(dif'*(D.*dif)) + Beta_reg/p * norm(Phi.times(prim), p)^p - Beta_reg/p * norm(Phi.times(x), p)^p;
    prim_sh = Phi.times(reshape(prim,xSz));
    
    Delta = dot(dif(:),grad_smooth(:)) + 1/(2*Alpha_k)*(dif'*(D.*dif)) + Beta_reg * SHpnorm(prim_sh, p) - Beta_reg * x_sh_pnorm;
    Psi = -fv;  % Psi is the objective function of the dual problem (ie, the opposite of the objective funtion fv that arise
                % from the min problem equivalent to max dual problem with Psi): don't mistake the optimization problem with the stopping criterion
    if  Delta <= eta*Psi || cont >= maxinnerit
        loop = false;
    end
    
    fprintf('\n         lambda = %e, Delta = %e, eta*Psi=%e\n', lam, Delta, eta*Psi)
    
    
end

cont = cont - 1;


return














