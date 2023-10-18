function N = cellnorm(x,p)
%CELLNORM Norm of a cell array
N = 0;
L = length(x);

if p == 1
    for l = 1:L
        N = N + norm(vec(x{l}),1);
    end
else
    for l = 1:L
        N = N + sum(abs(x{l}).^p,"all");
    end
    N = N^(1/p);
end
end

