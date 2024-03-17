function [zero, x, k, ae, re] = broydensys( ...
    fun, B0, x0, kmax, tolerance ...
)
    % init env
    x = {x0};
    fx = {fun(x0)};
    dx = {0};
    B = {B0};
    ae = Inf;
    re = Inf;

    k = 1;
    while re(k) > tolerance && k < kmax
        dx{k + 1} = - B{k} \ fx{k};
        x{k + 1} = x{k} + dx{k + 1};

        fx{k + 1} = fun(x{k + 1});

        df = fx{k + 1} - fx{k};
        B{k + 1} = B{k} + ((df - B{k} * dx{k + 1}) * dx{k + 1}') / (dx{k + 1}' * dx{k + 1});

        ae(k + 1) = norm(dx{k + 1});
        re(k + 1) = norm(dx{k + 1}) / norm(x{k});

        k = k + 1;
    end


    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    zero = x{k - 1};
end
