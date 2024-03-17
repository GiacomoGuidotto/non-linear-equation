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

%% Find zeros
% termination criteria
kmax = 1e3;
tolerance = 1e-8;

% initial points
x0A = [1; 1];
x0B = [-1; 1];

%% Newton method for both

[zeroA, xA, kA, dxA, aeA, reA] = newtonsys( ...
    @sys, @jac, x0A, kmax, tolerance ...
);

[zeroB, xB, kB, dxB, aeB, reB] = newtonsys( ...
    @sys, @jac, x0B, kmax, tolerance ...
);

%% Plot Newton method

f = figure();
f.Name = 'Non-linear system: Newton-Raphson method';
f.NumberTitle = 'off';
f.Position = [0, 0, 1000, 1500];

% Plot functions

subplot(3, 1, 1);
surf(X1, X2, Y1, 'FaceColor', 'r', 'EdgeColor', 'none');
hold on;
surf(X1, X2, Y2, 'FaceColor', 'b', 'EdgeColor', 'none');
hold on;
mesh(X1, X2, Z0);
title('3D plot of the system');
xlabel('x1'); ylabel('x2'); zlabel('y');

% Draw x lines and the zero points, green for the first point, black for the second

z_values = linspace(-10, 10, 2);
for i = 1:kA
    x_values = xA{i}(1) * ones(size(z_values));
    y_values = xA{i}(2) * ones(size(z_values));

    plot3(x_values, y_values, z_values, 'g', 'LineWidth', 2);
end
scatter3(zeroA(1), zeroA(2), 0, 'g', 'filled');

for i = 1:kB
    x_values = xB{i}(1) * ones(size(z_values));
    y_values = xB{i}(2) * ones(size(z_values));

    plot3(x_values, y_values, z_values, 'k', 'LineWidth', 2);
end
scatter3(zeroB(1), zeroB(2), 0, 'k', 'filled');

% Plot the two residual

subplot(3, 1, 2);
semilogy(1:kA, reA, '.-', 1:kB, reB, '.-');
title("residual decay", 'interpreter', 'latex');
xlabel("k"); ylabel("residual");
yline(tolerance, ':', 'tolerance');
legend(['$x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ')^T$'], ...
 ['$x0 = (', num2str(x0B(1)), ', ', num2str(x0B(1)), ')^T$'], ...
 'Location', 'northeast', 'interpreter', 'latex')

% Compute the convergence ratios

ratio1 = reA(2:end) ./ reA(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reA(2:end) ./ (reA(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

subplot(3, 1, 3);
semilogy(1:length(ratio1), ratio1, '.-', ...
        1:length(ratio2), ratio2, '.-');
title('ratios convergence (for the first point)', 'interpreter', 'latex');
legend('$\frac{d_k}{d_{k-1}}$', ...
        '$\frac{d_k}{d_{k-1}^2}$', ...
         'interpreter', 'latex');
xlabel("k"); ylabel("ratio");


% Log the final values

disp("Non-linear systems: Netwon-Raphson method");
disp(['first zero = (', num2str(zeroA(1)), ', ', num2str(zeroA(2)), ') with ', num2str(kA), ' iterations']);
disp(['second zero = (', num2str(zeroB(1)), ', ', num2str(zeroB(2)), ') with ', num2str(kB), ' iterations']);

%% Broyden method for both

B0 = eye(2);

[zero1, ~, k1] = broydensys( ...
    @sys, B0, x0A, kmax, tolerance ...
);

[zero2, ~, k2] = broydensys( ...
    @sys, B0, x0B, kmax, tolerance ...
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

