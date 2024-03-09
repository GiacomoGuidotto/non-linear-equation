function [final_guess, guesses, interations, residual] = newton( ...
    fun, fun_prime, x0, kmax, tolerance ...
)
    newton_function = @(x) x - fun(x) / fun_prime(x);

    [final_guess, guesses, interations, residual] = picard( ...
        newton_function, x0, kmax, tolerance ...
    );

    k = interations;

    figure('Name','Newton method','NumberTitle','off');
    subplot(2, 1, 1);
    semilogy(1:k, guesses(1:k), '.-');
    title("approximation of p_m");
    xlabel("k"); ylabel("guess");

    subplot(2, 1, 2);
    semilogy(1:k, residual(1:k), '.-');
    title("residual decay");
    xlabel("k"); ylabel("residual");
end