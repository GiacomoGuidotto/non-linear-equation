function [zero, x, k, re] = picard( ...
    fun, x0, kmax, tolerance ...
)
    % init env
    x = x0;
    re = abs(x0 - fun(x0)) / x0;

    % Picard iteration
    k = 1;
    while re(k) > tolerance && k < kmax
        x(k + 1) = fun(x(k));

        re(k + 1) = abs(x(k + 1) - x(k)) / x(k);

        k = k + 1;
    end

    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x(k - 1);
end
