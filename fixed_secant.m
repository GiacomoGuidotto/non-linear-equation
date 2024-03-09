function [final_guess, guesses, interations, residual] = fixed_secant( ...
    fun, x0, x1, kmax, tolerance ...
)
    slope = (fun(x1) - fun(x0)) / (x1 - x0);
    fixed_secant_function = @(x) x - fun(x) / slope;

    [final_guess, guesses, interations, residual] = picard( ...
        fixed_secant_function, x0, kmax, tolerance ...
    );

    k = interations;

    figure('Name','Fixed Secant method','NumberTitle','off');
    subplot(2, 1, 1);
    semilogy(1:k, guesses(1:k), '.-');
    title("approximation of p_m");
    xlabel("k"); ylabel("guess");

    subplot(2, 1, 2);
    semilogy(1:k, residual(1:k), '.-');
    title("residual decay");
    xlabel("k"); ylabel("residual");
end