function [final_guess, guesses, k, residual] = picard( ...
    fun, x0, kmax, tolerance ...
)
    % x0 is already root
    if fun(x0) == 0
        final_guess = x0; k = 1; guesses = x0; residual = 0;
        return;
    end

    % init env
    k = 1;
    guesses = zeros(kmax, 1);
    residual = zeros(kmax, 1);

    % initial values
    guesses(k) = x0;
    residual(k) = abs(x0 - fun(x0)) / x0;

    % Picard iteration
    while k < kmax && residual(k) > tolerance
        new_guess = fun(guesses(k));

        guesses(k + 1) = new_guess;
        residual(k + 1) = abs(new_guess - guesses(k)) / abs(new_guess);

        k = k + 1;
    end

    % trim vectors
    guesses = guesses(1:k);
    residual = residual(1:k);

    % termination criteria check
    if k == kmax
        disp('Warning: Maximum iterations reached');
    end

    % return
    final_guess = guesses(k - 1);
end
