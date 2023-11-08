function C = SH2feas(C, bound)
    % Constrain (shearlet) coefficients C to be inside [-bound, bound]
    for L = 1:length(C)
        b = max(C{L}, -bound);
        C{L} = min(b, bound);
    end
end
