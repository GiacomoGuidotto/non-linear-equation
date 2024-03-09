function [final_guess, guesses, interations, residual] = secant( ...
    fun, x0, x1, kmax, tolerance ...
)
    prevx = x0;
    function new_x = secant_function(x)
        slope = (fun(x) - fun(prevx)) / (x - prevx);
        new_x = x - fun(x) / slope;
        prevx = x;
    end

    [final_guess, guesses, interations, residual] = picard( ...
        @secant_function, x1, kmax, tolerance ...
    );

    k = interations;

    figure('Name','Secant method','NumberTitle','off');
    subplot(2, 1, 1);
    semilogy(1:k, guesses(1:k), '.-');
    title("approximation of p_m");
    xlabel("k"); ylabel("guess");

    subplot(2, 1, 2);
    semilogy(1:k, residual(1:k), '.-');
    title("residual decay");
    xlabel("k"); ylabel("residual");
end