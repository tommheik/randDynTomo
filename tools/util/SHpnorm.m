function pnorm = SHpnorm(coeff,p)
% SHPNORM is the l^p norm of the cylindrical shearlet coefficients (or any
% other coefficients stored in cell arrays)
% Input
%   coeff       Shearlet coefficients stored in a cell array for each level
%   p           Exponent p
% 
% Output
%   pnorm      1/p * \| SH(f) \|_p^p
% Note: We do not take the p'th root because it is very slow! Therefore we
% shouldn't use the norm-function either.
%
% T. Heikkilä    2022

S = zeros(1,length(coeff));
for k = 1:length(coeff)
    S(k) = sum(abs(coeff{k}).^p,'all');
end
pnorm = sum(S) / p;

