function [final_guess, guesses, k, residual] = picard( ...
    fun, x0, kmax, tolerance ...
)
    % Check if root was given
    if fun(x0) == 0
        final_guess = x0;
        k = 1;
        guesses = x0;
        residual = 0;
        return;
    end


    guesses = zeros(kmax, 1);
    residual = zeros(kmax, 1);

    guesses(1) = x0;
    residual(1) = abs(x0 - fun(x0)) / x0;

    % Perform Picard iteration
    k = 1;
    while k < kmax && residual(k) > tolerance
        new_guess = fun(guesses(k));

        guesses(k + 1) = new_guess;
        residual(k + 1) = abs(new_guess - guesses(k)) / abs(new_guess);

        k = k + 1;
    end

    % Check if maximum iterations reached
    if k == kmax
        disp('Warning: Maximum iterations reached.');
    end

    disp(['k', num2str(k)]);
    % Return final guess
    final_guess = guesses(k - 1);
end
