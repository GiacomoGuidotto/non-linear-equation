% Clear environment

clear
close all
clc

%% Part 2 - Non linear systems

%% System definition


%% Find zeros
% termination criteria
kmax = 1e3;
tolerance = 1e-8;

% initial points
x0 = [1, 1]';
x1 = [-1, 1]';

% Newton method for both

[zero1, ~, k1] = newtonsys( ...
    @sys, @jac, x0, kmax, tolerance ...
);

[zero2, ~, k2] = newtonsys( ...
    @sys, @jac, x1, kmax, tolerance ...
);

disp("Netwon method");
disp(['zero1 = [', num2str(zero1(1)), ' ', num2str(zero1(2)), ']']);
disp(['k1 = ', num2str(k1)]);
disp(['zero2 = [', num2str(zero2(1)), ' ', num2str(zero2(2)), ']']);
disp(['k2 = ', num2str(k2)]);

% Broyden method for both

B0 = eye(2);

[zero1, ~, k1] = broydensys( ...
    @sys, B0, x0, kmax, tolerance ...
);

[zero2, ~, k2] = broydensys( ...
    @sys, B0, x1, kmax, tolerance ...
);

disp("Broyden method");
disp(['zero1 = [', num2str(zero1(1)), ' ', num2str(zero1(2)), ']']);
disp(['k1 = ', num2str(k1)]);
disp(['zero2 = [', num2str(zero2(1)), ' ', num2str(zero2(2)), ']']);
disp(['k2 = ', num2str(k2)]);



%% Appendix
% Tested system

function y = sys(x)
    y1 = x(1)^2 + x(2)^2 - 1;
    y2 = sin(pi/2 * x(1)) + x(2)^3;
    y = [y1; y2];
end

function J = jac(x)
    J(1, 1) = 2 * x(1);
    J(1, 2) = 2 * x(2);
    J(2, 1) = (pi/2) * cos(pi/2 * x(1));
    J(2, 2) = 3 * x(2)^2;
end

