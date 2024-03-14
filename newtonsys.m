function [zero, x, k, re] = newtonsys( ...
    fun, jac, x0, kmax, tolerance ...
)
    % init env
    x = {x0};
    re = Inf;

    k = 1;
    while re(k) > tolerance && k < kmax
        dx = - jac(x{k}) \ fun(x{k});
        x{k + 1} = x{k} + dx;

        re(k + 1) = norm(dx);

        k = k + 1;
    end


    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x{k - 1};
end
