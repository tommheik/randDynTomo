function bregDistance = bregL1(Mu, Mv)
%BREGL1 Symmetric Bregman distance for p=1
%   Mu and Mv should be the input vectors u and v already transformed using
%   the representation system M
%   This implementation matches the choice for the subdifferential which is
%   r_f = sgn(Mf), i.e. 0 when Mf(k) = 0.
%
% Note: if one of the vectors is 0, then l1-norm should be used!
%
% T. Heikkilä   2023

lvl = length(Mu);
if length(Mv) ~= lvl
    error('Decomposition levels should agree! Now L_u = %d but L_v = %d', lvl, length(Mv))
end

bregDistance = 0;
for l = 1:lvl
    u = Mu{l};
    v = Mv{l};
    r = sign(u) - sign(v);
    rBool = (r(:) ~= 0);
    bregDistance = bregDistance + dot(r(rBool), u(rBool) - v(rBool));
end

