clc; clear; close all;

n = 300; % number of discretization points
s = linspace(0, 1, n)';

% Define 10 different initial profiles (dims)
dims = [
    exp(-5*s) + 1, ...
    sqrt(s) + 0.5*cos(3*pi*s), ...
    s.^2 + 0.2*sin(6*pi*s), ...
    mod(5*s, 1), ...
    log(s + 1e-2) + 3, ...
    sin(2*s+1) + 5, ...
    exp(-(s - 0.3).^2 / 0.01) + 0.5 * exp(-(s - 0.7).^2 / 0.02), ...
    cos(10*pi*s) + 2, ...
    double(s > 0.3 & s < 0.7), ...
    1 - abs(2*s - 1)
];

algNames = {'Ex2Alg3p3', 'Ex2CTC', 'Ex2DCZR', 'Ex2M', 'Ex2TVC', 'Ex2TQY', 'Ex2T'};

% Initialize results
for a = 1:length(algNames)
    for d = 1:size(dims, 2)
        results.(algNames{a})(d).iter = [];
        results.(algNames{a})(d).time = [];
    end
end

% Run experiments
for a = 1:length(algNames)
    alg = str2func(algNames{a});
    for d = 1:size(dims, 2)
        x0 = dims(:, d); % initial function profile
        tic;
        [~, k] = alg(x0);
        elapsed = toc;

        results.(algNames{a})(d).iter = k;
        results.(algNames{a})(d).time = elapsed;
    end
end

% ===== Generate LaTeX table =====
fprintf('\\begin{tabular}{l');
fprintf(' c', 2*size(dims, 2));
fprintf('}\n');
fprintf('\\hline\n');
fprintf('Algorithm ');
for d = 1:size(dims, 2)
    fprintf('& \\multicolumn{2}{c}{$n=%d$} ', n);
end
fprintf('\\\\\n');
fprintf(' ');
for d = 1:size(dims, 2)
    fprintf('& Iter & Time ');
end
fprintf('\\\\ \\hline\n');

for a = 1:length(algNames)
    fprintf('%s ', algNames{a});
    for d = 1:size(dims, 2)
        fprintf('& %d & %.4f ', results.(algNames{a})(d).iter, results.(algNames{a})(d).time);
    end
    fprintf('\\\\\n');
end

fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
