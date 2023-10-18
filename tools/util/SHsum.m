function W = SHsum(U, V, a, b)
%SHSUM W = aU + bV
    for l = 1:length(U)
        W{l} = a*U{l} + b*V{l};
    end
end

