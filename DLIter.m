clc; clear; close all;

% ------------------------------------------
% Iteration data for algorithms from Table T1
% Columns correspond to k = [10, 50, 100, 200, 300, 500, 700, 900, 1200, 1500]
% ------------------------------------------
perf_data = [
    % k=10   k=50   k=100   k=200   k=300   k=500   k=700   k=900   k=1200   k=1500
     60,     60,     60,     60,     60,     61,     61,     61,     61,     61;     % PASY Alg.
    911,   1705,   2200,   3001,   3001,   3001,   3001,   3001,   3001,   3001;    % CTC Alg.
    138,    146,    147,    148,    149,    149,    149,    149,    136,    136;    % DCZR Alg.
    110,    111,    113,    116,    118,    120,    123,    124,    126,    128;    % M Alg.
    427,   1092,   1530,   2140,   2589,   3001,   3001,   3001,   3001,   3001;    % TQY Alg.
    126,    236,    299,    373,    423,    494,    546,    588,    639,    681;    % TVC Alg.
    229,    236,    243,    252,    256,    262,    266,    269,    272,    274     % T Alg.
];

alg_names = {
    'PASY Alg.', 'CTC Alg.', 'DCZR Alg.', ...
    'M Alg.', 'TQY Alg.', 'TVC Alg.', 'T Alg.'
};

[n_alg, n_prob] = size(perf_data);

% ------------------------------------------
% Compute performance ratios
% ------------------------------------------
r = zeros(n_alg, n_prob);
for j = 1:n_prob
    best = min(perf_data(:, j));
    r(:, j) = perf_data(:, j) / best;
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
% Plot performance profile
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
title('Dolan and Moré Performance Profile (Iterations)', 'FontSize', 14);
legend('Location', 'southeast');
xlim([1, max(tau)]);
ylim([0, 1.05]);
set(gca, 'FontSize', 12);
