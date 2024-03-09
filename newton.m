function [final_guess, guesses, interations, residual] = newton( ...
    fun, fun_prime, x0, kmax, tolerance ...
)
    newton_function = @(x) x - fun(x) / fun_prime(x);

    [final_guess, guesses, interations, residual] = picard( ...
        newton_function, x0, kmax, tolerance ...
    );

    k = interations;

    f = figure();
    f.Name = 'Newton-Raphson method';
    f.NumberTitle = 'off';
    f.Position = [0, 0, 800, 900]

    subplot(3, 2, 1);
    semilogy(1:k, guesses, '.-');
    title("approximation of $p_m$", 'interpreter', 'latex');
    xlabel("k"); ylabel("guess");
    yline(final_guess, ':', 'final guess');

    subplot(3, 2, 3);
    semilogy(1:k, residual, '.-');
    title("residual decay", 'interpreter', 'latex');
    xlabel("k"); ylabel("residual");
    yline(tolerance, ':', 'tolerance');

    subplot(3, 2, 2);
    ratio_1 = residual(2:end) ./ residual(1:end-1);
    semilogy(2:k, ratio_1, '.-');
    title('ratio 1: $\frac{d_k}{d_{k-1}}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");

    subplot(3, 2, 4);
    ratio_2 = residual(2:end) ./ (residual(1:end-1) .^ 2)
    semilogy(2:k, ratio_2, '.-');
    title('ratio 2: $\frac{d_k}{d_{k-1}^2}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");

    subplot(3, 2, 6);
    golden = (1 + sqrt(5)) / 2
    ratio_3 = residual(2:end) ./ (residual(1:end-1) .^ golden)
    semilogy(2:k, ratio_3, '.-');
    title('ratio 3: $\frac{d_k}{d_{k-1}^{1.618}}$', 'interpreter', 'latex');
    xlabel("k"); ylabel("ratio");
end