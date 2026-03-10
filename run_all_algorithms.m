clc; clear; close all;

dims = [10, 50, 100, 200, 300, 500, 700, 900, 1200, 1500];
algNames = {'Ex1Alg3p3', 'Ex1CTC', 'Ex1DCZR', 'Ex1M', 'Ex1TVC', 'Ex1TQY', 'Ex1T'};

% Initialize results struct
for a = 1:length(algNames)
    for d = 1:length(dims)
        results.(algNames{a})(d).iter = [];
        results.(algNames{a})(d).time = [];
    end
end

% Run experiments
for a = 1:length(algNames)
    alg = str2func(algNames{a});
    for d = 1:length(dims)
        n = dims(d);
        tic;
        [~, k] = alg(n);
        elapsed = toc;
        
        results.(algNames{a})(d).iter = k;
        results.(algNames{a})(d).time = elapsed;
    end
end

% ===== Generate LaTeX table =====
fprintf('\\begin{tabular}{l');
fprintf(' c', 2*length(dims)); % 2 columns (iter, time) per dimension
fprintf('}\n');
fprintf('\\hline\n');
fprintf('Algorithm ');
for d = 1:length(dims)
    fprintf('& \\multicolumn{2}{c}{$n=%d$} ', dims(d));
end
fprintf('\\\\\n');
fprintf(' ');
for d = 1:length(dims)
    fprintf('& Iter & Time ');
end
fprintf('\\\\ \\hline\n');

for a = 1:length(algNames)
    fprintf('%s ', algNames{a});
    for d = 1:length(dims)
        fprintf('& %d & %.4f ', results.(algNames{a})(d).iter, results.(algNames{a})(d).time);
    end
    fprintf('\\\\\n');
end

fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
