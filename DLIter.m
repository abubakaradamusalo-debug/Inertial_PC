clc; clear; close all;

% ------------------------------------------
% Iteration data for algorithms from Table T1
% Columns correspond to k = [10, 50, 100, 200, 300, 500, 700, 900, 1200, 1500]
% ------------------------------------------
perf_data = [
    % Ex2Alg3p3
    60 56 71 68 74 81 72 73 61 68;
    % Ex2CTC
    453 437 417 426 475 509 414 473 430 426;
    % Ex2DCZR
    788 772 755 757 823 867 751 819 779 757;
    % Ex2M
    103 122 132 126 103 58 133 102 126 126;
    % Ex2TVC
    519 517 346 3001 338 242 1038 290 252 282;
    % Ex2TQY
    509 662 669 601 704 216 668 649 807 601;
    % Ex2T
    112 110 109 110 111 112 109 111 107 110
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
