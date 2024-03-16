function [zero, x, k, dx, re] = picard( ...
    kernel, x0, kmax, tolerance ...
)
    % init env
    x = x0;
    re = Inf(1, length(x));
    dx = zeros(1, length(x));

    % init vectors
    for k = 2:length(x)
        dx(k) = abs(x(k) - x(k - 1));
        re(k) = dx(k) / x(k - 1);
    end

    % Picard iteration
    k = length(x);
    while re(k) > tolerance && k < kmax
        x(k + 1) = kernel(x, k);

        dx(k + 1) = abs(x(k + 1) - x(k));
        re(k + 1) = dx(k + 1) / x(k);

        k = k + 1;
    end

    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x(k - 1);
end
