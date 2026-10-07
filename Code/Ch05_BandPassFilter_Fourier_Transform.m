%----------------------------------------------------------
%  Band-pass filter with Fourier transform computed by trapz
%  Optimized version with better visualization and efficiency
%----------------------------------------------------------
clc; clear; close all;

%% 1. Create a synthetic test signal
Fs   = 1000;                % sampling frequency (Hz)
T    = 1;                   % signal duration (s)
t    = (0:1/Fs:T-1/Fs)';    % time vector (column)
N    = numel(t);

% Signal components
f1 =  20;  A1 = 0.8;        % low-frequency component (to be rejected)
f2 = 120;  A2 = 0.5;        % mid-frequency component (wanted)
f3 = 250;  A3 = 0.4;        % high-frequency component (to be rejected)
noise = 0.3*randn(N,1);     % additive white noise

x = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t) + A3*sin(2*pi*f3*t) + noise;

%% 2. Angular-frequency axis
domega = 2*pi*0.25;         % angular-frequency resolution (rad/s)
omegaMax = 2*pi*(Fs/2);     % Nyquist angular frequency
omega = (-omegaMax:domega:omegaMax)';  % column vector (rad/s)
M = numel(omega);

%% 3. Compute Fourier transform X(omega) = \int x(t) e^(-j omega t) dt
dt = 1/Fs;                  % sampling period
X = zeros(M, 1, 'like', 1+1i);

% Precompute for efficiency
for k = 1:M
    phasor = exp(-1i*omega(k)*t);
    X(k) = trapz(t, x .* phasor);
end
% No scaling by dt needed here since trapz with time vector handles this

%% 4. Create band-pass filter
f_cl = 80;                  % lower cut-off frequency (Hz)
f_cu = 160;                 % upper cut-off frequency (Hz)
omega_cl = 2*pi*f_cl;       % lower angular cut-off (rad/s)
omega_cu = 2*pi*f_cu;       % upper angular cut-off (rad/s)

H = double((abs(omega) >= omega_cl) & (abs(omega) <= omega_cu));

% Visualise the filter
figure;
plot(omega/(2*pi), H, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Gain');
title('Ideal Band-Pass Filter Frequency Response');
grid on;
xlim([0, 300]);  % Show only positive frequencies for clarity

%% 5. Apply the filter
Y = X .* H;

%% 6. Inverse Fourier transform x(t) = (1/(2pi)) \int X(omega) e^(j omega t) domega
y = zeros(N, 1, 'like', 1+1i);

for n = 1:N
    phasor = exp(1i*omega*t(n));
    y(n) = trapz(omega, Y .* phasor);
end
y = real(y) / (2*pi);       % Scale by 1/(2pi) and take real part

%% 7. Plot results with improved visualization
figure;

% Time domain
subplot(3,1,1);
plot(t, x, 'b', 'DisplayName', 'Original'); hold on;
plot(t, y, 'r', 'LineWidth', 1.5, 'DisplayName', 'Filtered');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Location', 'best'); grid on;
title('Time-Domain Signals');
xlim([0, 0.2]);  % Show first 0.2 seconds for clarity

% Frequency domain (magnitude)
subplot(3,1,2);
f_axis = omega/(2*pi);
plot(f_axis, abs(X), 'b', 'LineWidth', 2.5,'DisplayName', 'Original spectrum'); hold on;
plot(f_axis, abs(Y), 'r', 'LineWidth', 1.5, 'DisplayName', 'Filtered spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude');
legend('Location', 'best'); grid on;
xlim([0, 300]);  % Show only positive frequencies
title('Magnitude Spectrum');

% Frequency domain (dB scale)
subplot(3,1,3);
db_X = 20*log10(abs(X)/max(abs(X)));
db_Y = 20*log10(abs(Y)/max(abs(Y)));
plot(f_axis, db_X, 'b', 'LineWidth', 2.5,'DisplayName', 'Original spectrum'); hold on;
plot(f_axis, db_Y, 'r', 'LineWidth', 1.5, 'DisplayName', 'Filtered spectrum');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Location', 'best'); grid on;
xlim([0, 300]);
ylim([-80, 0]);  % Typical dB range
title('Normalized Spectrum (dB scale)');

%% 8. Performance analysis
fprintf('Filter Performance Analysis:\n');
fprintf('Original signal power: %.4f\n', var(x));
fprintf('Filtered signal power: %.4f\n', var(y));
fprintf('Noise reduction: %.2f dB\n', 10*log10(var(x)/var(noise)));

% Check if desired component is preserved
t_analysis = t(1:min(1000, N));  % Use first second for frequency analysis
[P_orig, f_orig] = pwelch(x(1:length(t_analysis)), [], [], [], Fs);
[P_filt, f_filt] = pwelch(y(1:length(t_analysis)), [], [], [], Fs);

desired_power_orig = trapz(f_orig((f_orig>=f_cl) & (f_orig<=f_cu)), ...
                           P_orig((f_orig>=f_cl) & (f_orig<=f_cu)));
desired_power_filt = trapz(f_filt((f_filt>=f_cl) & (f_filt<=f_cu)), ...
                           P_filt((f_filt>=f_cl) & (f_filt<=f_cu)));

fprintf('Power in passband (80-160 Hz):\n');
fprintf('  Original: %.4f\n', desired_power_orig);
fprintf('  Filtered: %.4f\n', desired_power_filt);
fprintf('  Preservation: %.1f%%\n', 100*desired_power_filt/desired_power_orig);