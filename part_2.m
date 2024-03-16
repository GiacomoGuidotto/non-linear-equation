% Clear environment

clear
close all
clc

%% Part 2 - Non linear systems

%% System definition
f1 = @(x, y) x^2 + y^2 - 1;
f2 = @(x, y) sin(pi/2 * x) + y^3;

% Define the range of x values
x1_range = -2:0.1:2;
x2_range = -2:0.1:2;
[X1, X2] = meshgrid(x1_range, x2_range);

% Calculate the corresponding y values
Y1 = zeros(size(X1));
Y2 = zeros(size(X1));
Z0 = zeros(size(X1));
for i = 1:numel(X1)
    Y1(i) = f1(X1(i), X2(i));
    Y2(i) = f2(X1(i), X2(i));
end

% Plot the functions
figure;
surf(X1, X2, Y1, 'FaceColor', 'r', 'EdgeColor', 'none');
hold on;
surf(X1, X2, Y2, 'FaceColor', 'b', 'EdgeColor', 'none');
hold on;
mesh(X1, X2, Z0);
title('3D plot of the system');
xlabel('x1'); ylabel('x2'); zlabel('y');

%% Find zeros
% termination criteria
kmax = 1e3;
tolerance = 1e-8;

% initial points
x0 = [1, 1]';
x1 = [-1, 1]';

% Newton method for both

[zero1, x, k1, dx, re] = newtonsys( ...
    @sys, @jac, x0, kmax, tolerance ...
);

[zero2, xalt, k2] = newtonsys( ...
    @sys, @jac, x1, kmax, tolerance ...
);

for i = 1:k1
    scatter3(x{i}(1), x{i}(2), 0, 'g', 'filled');
end
for i = 1:k2
    scatter3(xalt{i}(1), xalt{i}(2), 0, 'k', 'filled');
end


scatter3(zero1(1), zero1(2), 0, 'g', 'filled');
scatter3(zero2(1), zero2(2), 0, 'k', 'filled');


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

