clc;
clear;
close all;

% Define number of grid points
N = 1001;

% Define spatial step size
h = 10 / (N - 1);

% Define maximum wave speed
u_max = 1;

% Define computational domain
xi = -4.995:0.01:4.995;  % Fine grid for numerical solution
x = linspace(-5, 5, 1000); % Coarse grid for plotting

% Length of the computational domain
a = length(xi);

% Initial condition
u0 = zeros(1, a);   % Initialize solution array with zeros
u0(x > -1 & x < 0) = -1;  % Define initial condition for -1 < x < 0
u0(x > 0 & x < 1) = 1;    % Define initial condition for 0 < x < 1

% Duplicate initial condition for different methods
u0(2, :) = u0(1, :);  % For method 2
u0(3, :) = u0(1, :);  % For method 3

% Define time step size based on stability condition (CFL condition)
t = 0.5 * h / u_max;

% Number of time steps to reach t = 2 seconds
nt = floor(2 / t);

% Initialize solution matrix for all time steps
u_p = u0;

% Define flux function
f = @(u) 0.5 * u^2;

% Define slope limiters
Mind = @(r) max(0, min(r, 1));    % Minmod limiter
sbee = @(r) max(max(0, min(2*r, 1)), min(r, 2));  % Superbee limiter
non = @(r) 1;  % No limiter (constant slope)

% Store limiters in a cell array
Lim = {non, Mind, sbee};

% Main loop for all limiters (methods)
for j = 1:3
    % Time-stepping loop
    for n = 1:nt
        % Copy current solution
        u1(j, :) = u_p(j, :);

        % Loop over interior grid points (excluding boundaries)
        for i = 3:a-2
            % Compute left and right states at cell interfaces
            uL_1 = U(u_p(j, i-1), x(i-1) + 0.5 * h, x(i-1), Lim{j}, r(u_p(j, i-2), u_p(j, i-1), u_p(j, i)), s(u_p(j, i), u_p(j, i-1), h));
            uR_1 = U(u_p(j, i), x(i-1) + 0.5 * h, x(i), Lim{j}, r(u_p(j, i-1), u_p(j, i), u_p(j, i+1)), s(u_p(j, i+1), u_p(j, i), h));
            uL_2 = U(u_p(j, i), x(i) + 0.5 * h, x(i), Lim{j}, r(u_p(j, i-1), u_p(j, i), u_p(j, i+1)), s(u_p(j, i+1), u_p(j, i), h));
            uR_2 = U(u_p(j, i+1), x(i) + 0.5 * h, x(i+1), Lim{j}, r(u_p(j, i), u_p(j, i+1), u_p(j, i+2)), s(u_p(j, i+2), u_p(j, i+1), h));

            % Compute numerical fluxes at interfaces
            flux_left = g_flux(uL_1, uR_1, f);
            flux_right = g_flux(uL_2, uR_2, f);

            % Update solution using Godunov's method
            u1(j, i) = u_p(j, i) - (t / h) * (flux_right - flux_left);
        end

        % Update solution for the next time step
        u_p(j, :) = u1(j, :);
    end
end

% Plot the result
plot(x, u1, 'LineWidth', 2);
title('Solution at t = 2 s using Godunov''s Method');
xlabel('x');
ylabel('u(x, t)');

% ---------------- Helper Functions ----------------

% Numerical flux function
function F = g_flux(uL, uR, f)
    % Compute numerical flux based on left and right states
    F = (uL > uR) * (max(f(uL), f(uR))) + ...
        (uL <= uR) * (min(f(uL), f(uR)));
end

% Slope ratio function
function a = r(uL, uM, uR)
    % Compute slope ratio between neighboring cells
    a = (uM - uL) / (uR - uM);
end

% Reconstruction function
function a = U(ui, x, xi, phi, r, s)
    % Reconstruct the solution using the slope limiter
    a = ui + phi(r) * s * (x - xi) / 2;
end

% Slope function
function a = s(uR, uL, h)
    % Compute slope between neighboring cells
    a = (uR - uL) / h;
end
