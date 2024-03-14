function zero = newton( ...
    fun, fun_prime, x0, kmax, tolerance ...
)
    newton_kernel = @(x, k) x(k) - fun(x(k)) / fun_prime(x(k));

    [zero, x, ~, residual] = picard( ...
        newton_kernel, x0, kmax, tolerance ...
    );

    ratio_1 = residual(2:end) ./ residual(1:end-1);

    ratio_2 = residual(2:end) ./ (residual(1:end-1) .^ 2);

    golden = (1 + sqrt(5)) / 2;
    ratio_3 = residual(2:end) ./ (residual(1:end-1) .^ golden);

    f = figure();
    f.Name = 'Newton-Raphson method';
    f.NumberTitle = 'off';
    f.Position = [0, 0, 800, 900];

    subplot(3, 2, 1);
    semilogy(x, '.-');
    title("approximation of $p_m$", 'interpreter', 'latex');
    xlabel("k"); ylabel("guess");
    yline(zero, ':', 'final guess');

    subplot(3, 2, 3);
    semilogy(residual, '.-');
    title("residual decay", 'interpreter', 'latex');
    xlabel("k"); ylabel("residual");
    yline(tolerance, ':', 'tolerance');

    subplot(3, 2, 2);
    semilogy(ratio_1, '.-');
    title('ratio 1: $\frac{d_k}{d_{k-1}}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");

    subplot(3, 2, 4);
    semilogy(ratio_2, '.-');
    title('ratio 2: $\frac{d_k}{d_{k-1}^2}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");

    subplot(3, 2, 6);
    semilogy(ratio_3, '.-');
    title('ratio 3: $\frac{d_k}{d_{k-1}^{1.618}}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");

    T = table(x, residual, ratio_1, ratio_2, ratio_3);
    disp(T);
end