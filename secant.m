function [zero, x, k, re] = secant( ...
    fun, x1, x2, kmax, tolerance ...
)
    % init env
    x = [x1, x2];
    y1 = fun(x1);
    re = abs(x2 - x1) / x1;

    % secant step
    k = 2;
    while re(k - 1) > tolerance && k < kmax
        y2 = fun(x(k));
        dx = - y2 * (x(k) - x(k - 1)) / (y2 - y1);
        x(k + 1) = x(k) + dx;

        re(k) = abs(x(k + 1) - x(k)) / x(k);

        k = k + 1;
        y1 = y2;
    end

    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x(k - 1);
end
