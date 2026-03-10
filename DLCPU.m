clc; clear; close all;

% ------------------------------------------
% CPU time data (in seconds) from Table T1
% Columns correspond to k = [10, 50, 100, 200, 300, 500, 700, 900, 1200, 1500]
% ------------------------------------------
cpu_time = [
    % Ex2Alg3p3
    22.7229 55.4162 71.0202 6.9553 32.2091 9.3411 79.3251 106.5810 2.9473 1.2339;
    % Ex2CTC
    2.7742 35.5811 81.9207 24.9354 7.5032 3.1208 79.5237 191.5230 13.5576 4.7088;
    % Ex2DCZR
    211.7009 537.9738 622.4271 76.3601 208.8139 52.7935 613.6063 1018.3340 46.2995 18.2937;
    % Ex2M
    0.9082 21.4718 56.7725 12.7537 2.5241 0.5750 59.8354 70.1339 6.8345 2.3626;
    % Ex2TVC
    260.9169 492.8708 403.7673 342.8572 185.1360 36.5118 1279.9710 484.5339 14.3560 6.5552;
    % Ex2TQY
    10.2309 87.8420 144.7966 34.1563 29.6606 2.0308 146.0535 302.0973 26.7527 8.0994;
    % Ex2T
    20.4725 59.1922 77.3744 11.2294 24.8260 7.5320 82.0636 119.0088 5.8723 2.2104
];

alg_names = {
    'PASY Alg.', 'CTC Alg.', 'DCZR Alg.', ...
    'M Alg.', 'TQY Alg.', 'TVC Alg.', 'T Alg.'
};

[n_alg, n_prob] = size(cpu_time);

% ------------------------------------------
% Compute performance ratios
% ------------------------------------------
r = zeros(n_alg, n_prob);
for j = 1:n_prob
    best = min(cpu_time(:, j));
    r(:, j) = cpu_time(:, j) / best;
end

% Tau values and cumulative distribution
tau = linspace(1, max(r(:)) + 0.1, 500);
rho = zeros(n_alg, length(tau));

for i = 1:n_alg
    for k = 1:length(tau)
        rho(i, k) = sum(r(i, :) <= tau(k)) / n_prob;
    end
end

% ------------------------------------------
% Colors and line styles
% ------------------------------------------
colors = lines(n_alg);
line_styles = {'-', '--', '-.', ':', '-', '--', '-.'};

% ------------------------------------------
% Plot performance profile for CPU time
% ------------------------------------------
figure('Units', 'pixels', 'Position', [100, 100, 750, 500], 'Color', 'w');
hold on; grid on;

for i = 1:n_alg
    plot(tau, rho(i, :), line_styles{i}, ...
        'Color', colors(i, :), 'LineWidth', 2.5, ...
        'DisplayName', alg_names{i});
end

xlabel('$\tau$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\rho_s(\tau)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Dolan and Moré Performance Profile (CPU Time)', 'FontSize', 14);
legend('Location', 'southeast');
xlim([1, max(tau)]);
ylim([0, 1.05]);
set(gca, 'FontSize', 12);
