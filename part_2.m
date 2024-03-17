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

%% Newton method

[zeroA, xA, kA, aeA, reA] = newtonsys( ...
    @sys, @jac, x0A, kmax, tolerance ...
);

[zeroB, xB, kB, aeB, reB] = newtonsys( ...
    @sys, @jac, x0B, kmax, tolerance ...
);

%% Plots

f = figure();
f.Name = 'Non-linear system: Newton-Raphson method';
f.NumberTitle = 'off';
f.Position = [0, 0, 800, 1300];

% Plot 1: approximations over the systems

subplot(3, 1, 1);
surf(X1, X2, Y1, 'FaceColor', 'r');
hold on;
surf(X1, X2, Y2, 'FaceColor', 'b');
hold on;
mesh(X1, X2, Z0);
title('zero approximations with Newton method', 'interpreter', 'latex');
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

% Plot 2: convergences errors

subplot(3, 1, 2);
semilogy(reA, '.-');
title(['residual decay with $x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ')^T$'], ...
    'interpreter', 'latex');
xlabel("k"); ylabel("residual");
yline(tolerance, ':', 'tolerance');

% Plot 3: convergence ratios

ratio1 = reA(2:end) ./ reA(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reA(2:end) ./ (reA(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

subplot(3, 1, 3);
semilogy(1:length(ratio1), ratio1, '.-', ...
        1:length(ratio2), ratio2, '.-');
title(['residual decay with $x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ')^T$'], ...
    'interpreter', 'latex');
legend('$\frac{d_k}{d_{k-1}}$', ...
        '$\frac{d_k}{d_{k-1}^2}$', ...
         'interpreter', 'latex');
xlabel("k"); ylabel("ratio");

% Print results in a table

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kA - 1))', aeA', reA', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'newtonsys_A.csv', 'Delimiter',';');

ratio1 = reB(2:end) ./ reB(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reB(2:end) ./ (reB(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kB - 1))', aeB', reB', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'newtonsys_B.csv', 'Delimiter',';');

% Log the final values

disp('zero approximations with Newton method');
disp(['first zero = (', num2str(zeroA(1)), ', ', num2str(zeroA(2)), ...
    ') from x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ') with ', num2str(kA), ' iterations']);
disp(['second zero = (', num2str(zeroB(1)), ', ', num2str(zeroB(2)), ...
    ') from x0 = (', num2str(x0B(1)), ', ', num2str(x0B(1)), ...
    ') with ', num2str(kB), ' iterations']);

%% Broyden method with B0 = I
B0 = eye(2);

[zeroA, xA, kA, aeA, reA] = broydensys( ...
    @sys, B0, x0A, kmax, tolerance ...
);

[zeroB, xB, kB, aeB, reB] = broydensys( ...
    @sys, B0, x0B, kmax, tolerance ...
);

%% Plots

f = figure();
f.Name = 'Non-linear system: Broyden method';
f.NumberTitle = 'off';
f.Position = [800, 0, 1600, 1300];

% Plot 1: approximations over the systems

subplot(3, 2, 1);
surf(X1, X2, Y1, 'FaceColor', 'r');
hold on;
surf(X1, X2, Y2, 'FaceColor', 'b');
hold on;
mesh(X1, X2, Z0);
title('zero approximations with Broyden method with $\textbf{B}^{(0)} = \textbf{I}$', ...
  'interpreter', 'latex');
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

% Plot 2: convergences errors

subplot(3, 2, 3);
semilogy(reA, '.-');
title(['residual decay with $x0 = (', ...
    num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ')^T$ and $\textbf{B}^{(0)} = \textbf{I}$' ...
    ], 'interpreter', 'latex');
xlabel("k"); ylabel("residual");
yline(tolerance, ':', 'tolerance');

% Plot 3: convergence ratios

ratio1 = reA(2:end) ./ reA(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reA(2:end) ./ (reA(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

subplot(3, 2, 5);
semilogy(1:length(ratio1), ratio1, '.-', ...
        1:length(ratio2), ratio2, '.-');
title(['residual decay with $x0 = (', ...
    num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ')^T$ and $\textbf{B}^{(0)} = \textbf{I}$' ...
    ], 'interpreter', 'latex');
legend('$\frac{d_k}{d_{k-1}}$', ...
        '$\frac{d_k}{d_{k-1}^2}$', ...
         'interpreter', 'latex');
xlabel("k"); ylabel("ratio");

% Print results in a table

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kA - 1))', aeA', reA', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'broydensys_B=I_A.csv', 'Delimiter',';');

ratio1 = reB(2:end) ./ reB(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reB(2:end) ./ (reB(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kB - 1))', aeB', reB', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'broydensys_B=I_B.csv', 'Delimiter',';');

% Log the final values

disp('zero approximations with Broyden method with B0 = I');
disp(['first zero = (', num2str(zeroA(1)), ', ', num2str(zeroA(2)), ...
    ') from x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ') with ', num2str(kA), ' iterations']);
disp(['second zero = (', num2str(zeroB(1)), ', ', num2str(zeroB(2)), ...
    ') from x0 = (', num2str(x0B(1)), ', ', num2str(x0B(1)), ...
    ') with ', num2str(kB), ' iterations']);

%% Broyden method with B0 = 2I
B0 = 2 * eye(2);

[zeroA, xA, kA, aeA, reA] = broydensys( ...
    @sys, B0, x0A, kmax, tolerance ...
);

[zeroB, xB, kB, aeB, reB] = broydensys( ...
    @sys, B0, x0B, kmax, tolerance ...
);

%% Plots

% Plot 1: approximations over the systems

subplot(3, 2, 2);
surf(X1, X2, Y1, 'FaceColor', 'r');
hold on;
surf(X1, X2, Y2, 'FaceColor', 'b');
hold on;
mesh(X1, X2, Z0);
title('zero approximations with Broyden method with $\textbf{B}^{(0)} = 2\textbf{I}$', ...
  'interpreter', 'latex');
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

% Plot 2: convergences errors

subplot(3, 2, 4);
semilogy(reA, '.-');
title(['residual decay with $x0 = (', ...
    num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ')^T$ and $\textbf{B}^{(0)} = 2\textbf{I}$' ...
    ], 'interpreter', 'latex');
xlabel("k"); ylabel("residual");
yline(tolerance, ':', 'tolerance');

% Plot 3: convergences ratios

ratio1 = reA(2:end) ./ reA(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reA(2:end) ./ (reA(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

subplot(3, 2, 6);
semilogy(1:length(ratio1), ratio1, '.-', ...
        1:length(ratio2), ratio2, '.-');
title(['residual decay with $x0 = (', ...
    num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ')^T$ and $\textbf{B}^{(0)} = 2\textbf{I}$' ...
    ], 'interpreter', 'latex');
legend('$\frac{d_k}{d_{k-1}}$', ...
        '$\frac{d_k}{d_{k-1}^2}$', ...
         'interpreter', 'latex');
xlabel("k"); ylabel("ratio");

% Print results in a table

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kA - 1))', aeA', reA', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'broydensys_B=2I_A.csv', 'Delimiter',';');

ratio1 = reB(2:end) ./ reB(1:end-1);
ratio1 = [Inf ratio1];

ratio2 = reB(2:end) ./ (reB(1:end-1) .^ 2);
ratio2 = [Inf ratio2];

varNames = ["k", "increment (normalized)", "relative increment (normilized)", "ratio p = 1", "ratio p = 2"];
T = table((0:(kB - 1))', aeB', reB', ratio1', ratio2', 'VariableNames', varNames);
writetable(T, 'broydensys_B=2I_B.csv', 'Delimiter',';');

% Log the final values

disp('zero approximations with Broyden method with B0 = 2I');
disp(['first zero = (', num2str(zeroA(1)), ', ', num2str(zeroA(2)), ...
    ') from x0 = (', num2str(x0A(1)), ', ', num2str(x0A(1)), ...
    ') with ', num2str(kA), ' iterations']);
disp(['second zero = (', num2str(zeroB(1)), ', ', num2str(zeroB(2)), ...
    ') from x0 = (', num2str(x0B(1)), ', ', num2str(x0B(1)), ...
    ') with ', num2str(kB), ' iterations']);

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
