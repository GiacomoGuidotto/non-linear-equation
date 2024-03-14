function [zero, x, k, re] = broydensys( ...
    fun, B0, x0, kmax, tolerance ...
)
    % init env
    x = {x0};
    B = {B0};
    re = Inf;

    k = 1;
    while re(k) > tolerance && k < kmax
        dx = - B{k} \ fun(x{k});
        x{k + 1} = x{k} + dx;

        df = fun(x{k + 1}) - fun(x{k});
        B{k + 1} = B{k} + ((df - B{k} * dx) * dx') / (dx' * dx);

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
