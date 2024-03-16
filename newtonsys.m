function [zero, x, k, dx, re] = newtonsys( ...
    fun, jac, x0, kmax, tolerance ...
)
    % init env
    x = {x0};
    dx = {0};
    re = Inf;

    k = 1;
    while re(k) > tolerance && k < kmax
        dx{k + 1} = - jac(x{k}) \ fun(x{k});
        x{k + 1} = x{k} + dx{k + 1};

        re(k + 1) = norm(dx{k + 1});

        k = k + 1;
    end


    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x{k - 1};
end
