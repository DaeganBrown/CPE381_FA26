% Numerical Fourier Transform of u(t) - u(t-2)
clear all; close all; clc;

% Parameters
T_max = 5;          % Maximum time (should cover the signal duration)
N_t = 10000;        % Number of time points
Omega_max = 20;     % Maximum frequency to compute
N_Omega = 1000;     % Number of frequency points

% Time vector
t = linspace(-T_max, T_max, N_t);
dt = t(2) - t(1);

% Define the function f(t) = u(t) - u(t-2)
f = zeros(size(t));
f(t >= 0 & t < 2) = 1;  % Pulse from t=0 to t=2

% Frequency vector
Omega = linspace(-Omega_max, Omega_max, N_Omega);

% Initialize Fourier transform
F_Omega = zeros(size(Omega));

% Numerical integration for each frequency
for k = 1:length(Omega)
    integrand = f .* exp(-1j * Omega(k) * t);
    F_Omega(k) = trapz(t, integrand);  % Trapezoidal integration
end

% Analytical solution for comparison
F_analytical = (1 - exp(-1j * 2 * Omega)) ./ (1j * Omega);
% Handle Omega = 0 case (limit = 2)
F_analytical(Omega == 0) = 2;

% Plot results
figure('Position', [100, 100, 1500, 800]);

% Plot 1: Time domain signal
subplot(2,3,1);
plot(t, f, 'b', 'LineWidth', 2);
xlabel('Time t');
ylabel('f(t)');
title('Time Domain: u(t) - u(t-2)');
xlim([-1, 3]);
ylim([-0.2, 1.2]);
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('Time','Interpreter', 'latex');
set(gca, 'FontSize', 15);
grid on;
grid minor;

% Plot 2: Magnitude of Fourier transform
subplot(2,3,2);
plot(Omega, abs(F_Omega), 'r', 'LineWidth', 2);
hold on;
plot(Omega, abs(F_analytical), 'k--', 'LineWidth', 1);
xlabel('Frequency \Omega', 'Interpreter','latex');
ylabel('|F(\Omega)|');
title('Magnitude Spectrum');
legend('Numerical', 'Analytical', 'Location', 'best');
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('$\Omega$','Interpreter', 'latex');
set(gca, 'FontSize', 15);
grid on;
grid minor;


% Plot 3: Phase of Fourier transform
subplot(2,3,3);
plot(Omega, angle(F_Omega), 'g', 'LineWidth', 2);
hold on;
plot(Omega, angle(F_analytical), 'k--', 'LineWidth', 1);
xlabel('Frequency \Omega');
ylabel('angle F(\Omega)');
title('Phase Spectrum');
legend('Numerical', 'Analytical', 'Location', 'best');
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('$\Omega$','Interpreter', 'latex');
set(gca, 'FontSize', 15);
grid on;
grid minor;

% Plot 4: Real part
subplot(2,3,4);
plot(Omega, real(F_Omega), 'm', 'LineWidth', 2);
hold on;
plot(Omega, real(F_analytical), 'k--', 'LineWidth', 1);
xlabel('Frequency \Omega');
ylabel('Re[F(\Omega)]');
title('Real Part');
legend('Numerical', 'Analytical', 'Location', 'best');
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('$\Omega$','Interpreter', 'latex');
set(gca, 'FontSize', 15);
grid on;
grid minor;

% Plot 5: Imaginary part
subplot(2,3,5);
plot(Omega, imag(F_Omega), 'c', 'LineWidth', 2);
hold on;
plot(Omega, imag(F_analytical), 'k--', 'LineWidth', 1);
xlabel('Frequency \Omega');
ylabel('Im[F(\Omega)]');
title('Imaginary Part');
legend('Numerical', 'Analytical', 'Location', 'best');
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('$\Omega$','Interpreter', 'latex');
set(gca, 'FontSize', 15);
grid on;
grid minor;

% Plot 6: Error comparison
subplot(2,3,6);
error = abs(F_Omega - F_analytical);
semilogy(Omega, error, 'k', 'LineWidth', 2);
xlabel('Frequency \Omega');
ylabel('Absolute Error');
title('Numerical vs Analytical Error');
set(gca, 'XColor', [0, 0, 0], 'YColor', [0, 0, 0], 'TickDir', 'out');
xaxis = get(gca, 'XAxis');
xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis');
yaxis.TickLabelInterpreter = 'latex';
xlabel('$\Omega$','Interpreter', 'latex');
set(gca, 'FontSize', 10);
grid on;
grid minor;

% Display numerical results at specific frequencies
fprintf('Fourier Transform Results:\n');
fprintf('Omega\t\tNumerical\t\tAnalytical\t\tError\n');
fprintf('------------------------------------------------\n');

test_Omegas = [0, 1, 2, 5, 10];
for i = 1:length(test_Omegas)
    [~, idx] = min(abs(Omega - test_Omegas(i)));
    fprintf('%.1f\t', Omega(idx));
    fprintf('%.6f%+.6fj\t', real(F_Omega(idx)), imag(F_Omega(idx)));
    fprintf('%.6f%+.6fj\t', real(F_analytical(idx)), imag(F_analytical(idx)));
    fprintf('%.2e\n', error(idx));
end

% Special check at Omega=0
[~, zero_idx] = min(abs(Omega));
fprintf('\nAt Omega=0:\n');
fprintf('Numerical: %.6f%+.6fj\n', real(F_Omega(zero_idx)), imag(F_Omega(zero_idx)));
fprintf('Analytical: %.6f%+.6fj\n', real(F_analytical(zero_idx)), imag(F_analytical(zero_idx)));
fprintf('Expected area: 2.000000\n');