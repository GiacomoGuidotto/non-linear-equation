function zero = secant( ...
    fun, x0, x1, kmax, tolerance ...
)
    y1 = fun(x0);

    function xnew = secant_kernel(x, k)
        y2 = fun(x(k));
        dx = - y2 * (x(k) - x(k - 1)) / (y2 - y1);
        xnew = x(k) + dx;
        y1 = y2;
    end


    [zero, x, ~, residual] = picard( ...
        @secant_kernel, [x0, x1], kmax, tolerance ...
    );

    ratio_1 = residual(2:end) ./ residual(1:end-1);
    ratio_1 = [Inf ratio_1];

    ratio_2 = residual(2:end) ./ (residual(1:end-1) .^ 2);
    ratio_2 = [Inf ratio_2];

    golden = (1 + sqrt(5)) / 2;
    ratio_3 = residual(2:end) ./ (residual(1:end-1) .^ golden);
    ratio_3 = [Inf ratio_3];

    f = figure();
    f.Name = 'Secant method';
    f.NumberTitle = 'off';
    f.Position = [1000, 0, 500, 900];

    subplot(3, 1, 1);
    semilogy(x, '.-');
    title("approximation of $p_m$", 'interpreter', 'latex');
    xlabel("k"); ylabel("guess");
    yline(zero, ':', 'final guess');

    subplot(3, 1, 2);
    semilogy(residual, '.-');
    title("residual decay", 'interpreter', 'latex');
    xlabel("k"); ylabel("residual");
    yline(tolerance, ':', 'tolerance');

    subplot(3, 1, 3);
    semilogy(1:length(ratio_1), ratio_1, '.-', ...
            1:length(ratio_2), ratio_2, '.-',...
            1:length(ratio_3), ratio_3, '.-');
    title('ratios convergence', 'interpreter', 'latex');
    legend('$\frac{d_k}{d_{k-1}}$', ...
            '$\frac{d_k}{d_{k-1}^2}$', ...
            '$\frac{d_k}{d_{k-1}^{1.618}}$', ...
             'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");
end