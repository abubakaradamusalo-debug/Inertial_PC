clc; clear; close all;

% ------------------------------------------
% CPU time data (in seconds) from Table T1
% Columns correspond to k = [10, 50, 100, 200, 300, 500, 700, 900, 1200, 1500]
% ------------------------------------------
cpu_time = [
    % k=10    k=50    k=100   k=200   k=300   k=500    k=700    k=900    k=1200   k=1500
    0.0043, 0.0058, 0.0127, 0.0559, 0.1183, 0.4781, 1.1066, 1.7807, 3.2445, 5.1478;   % PASY Alg.
    0.0280, 0.0806, 0.2344, 1.2057, 2.5691, 10.4383, 24.2785, 39.2131, 86.2253, 141.4386; % CTC Alg.
    0.0104, 0.0159, 0.0266, 0.1203, 0.2559, 0.8410, 1.7206, 2.6949, 4.2965, 6.3390;    % DCZR Alg.
    0.0067, 0.0198, 0.0220, 0.1194, 0.2214, 0.6813, 1.4312, 2.2684, 4.0047, 6.0685;    % M Alg.
    0.0355, 0.1862, 0.5769, 3.7975, 9.5319, 37.5968, 75.9712, 120.0453, 212.4845, 307.3920; % TQY Alg.
    0.0096, 0.0274, 0.0615, 0.3014, 0.6948, 2.8565, 6.4855, 10.6442, 19.9237, 31.5989;  % TVC Alg.
    0.0109, 0.0211, 0.0460, 0.2081, 0.4350, 1.4522, 3.1465, 5.0738, 8.6756, 12.5596      % T Alg.
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
