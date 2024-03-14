function zero = newton( ...
    fun, fun_prime, x0, kmax, tolerance ...
)
    newton_kernel = @(x, k) x(k) - fun(x(k)) / fun_prime(x(k));

    [zero, x, ~, residual] = picard( ...
        newton_kernel, x0, kmax, tolerance ...
    );

    ratio_1 = residual(2:end) ./ residual(1:end-1);
    ratio_1 = [Inf ratio_1];

    ratio_2 = residual(2:end) ./ (residual(1:end-1) .^ 2);
    ratio_2 = [Inf ratio_2];

    golden = (1 + sqrt(5)) / 2;
    ratio_3 = residual(2:end) ./ (residual(1:end-1) .^ golden);
    ratio_3 = [Inf ratio_3];

    f = figure();
    f.Name = 'Newton-Raphson method';
    f.NumberTitle = 'off';
    f.Position = [0, 0, 500, 900];

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

    varNames = ["x", "residual", "ratio p = 1", "ratio p = 2", "ratio p = ϕ"];
    T = table(x', residual', ratio_1', ratio_2', ratio_3', 'VariableNames', varNames);
    writetable(T, 'newton.csv', 'Delimiter',';');
end