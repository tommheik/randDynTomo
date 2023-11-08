function W = SHsum(U, V, a, b)
    % SHSUM W = aU + bV
    % Linear combination of cell arrays
    Llen = length(U);
    W  = cell(1,Llen);
    for L = 1:Llen
        W{L} = a*U{L} + b*V{L};
    end
end

