% Clear environment

clear
close all
clc

%% Part 2 - Non linear systems

%% System definition

% Interval to plot
x1 = linspace(-1, 1, 10);
x2 = linspace(-1, 1, 10);
[X1, X2] = meshgrid(x1, x2);

% Compute the corresponding values of y1 and y2
Y1 = X1.^2 + X2.^2 - 1;
Y2 = sin(pi/2 * X1) + X2.^3;

% Plot the surfaces
figure;
surf(X1, X2, Y1, 'FaceAlpha', 0.5);
hold on;
surf(X1, X2, Y2, 'FaceAlpha', 0.5);
xlabel('x1'); ylabel('x2'); zlabel('y');
title('3D plot of the system of equations');
legend('x1^2 + x2^2 - 1 = 0', 'sin(pi/2 * x1) + x2^3 = 0');
view(-40, 30);
grid on;
colorbar;
hold off;

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

