clc;
clear;
close all;

% Set random seed
seed = 97006855;
ss = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(ss);

% Parameters
ITR_max = 2000;
lambda  = 0.5;

para1.Itr_max = ITR_max; para1.gamma = 1; para1.sigma = 1;     para1.tau = 0.5; para1.rho = 1.5;
para2.Itr_max = ITR_max; para2.gamma = 1; para2.sigma = 0.001; para2.tau = 0.5; para2.rho = 1.8;

% Datasets
dataset = {'a1a.t','a2a.t','a3a.t','a4a.t','a5a.t','a6a.t','a7a.t','a8a.t','a9a.t'};

fid_tex = fopen('Rel.txt','w');

accuracy_results       = [];   % (#datasets) x 7
training_time_results  = [];   % (#datasets) x 7

for i = 1:length(dataset)
    [b, A] = libsvmread(dataset{i});
    [m, n] = size(A);
    fprintf('name=%s, m=%d, n=%d, lambda=%.2f\n', dataset{i}, m, n, lambda);

    % logistic regression loss handle (assumed available)
    fun = @(x) lr_loss(A, b, m, x, lambda);

    progress_r = [];  % rows = repeats, cols = [NI1 T1 G1 NI2 T2 G2 ... NI7 T7 G7]  (21 cols)
    acc_r      = [];  % rows = repeats, cols = [acc1 ... acc7]
    time_r     = [];  % rows = repeats, cols = [T1 ... T7]

    for repeat = 1:2
        % Initial point
        x0 = 90 * (ones(n,1));

        % Run 7 algorithms (T = training time)
        [T1, NI1, G1, x1] = IDFPM(fun, 'PASYAlg3p3', 3,  4, x0, para1);
        [T2, NI2, G2, x2] = IDFPM(fun, 'CTCAlg3p1', 3,  4, x0, para1);
        [T3, NI3, G3, x3] = IDFPM(fun, 'DCZRAlg3p1',   3,  6, x0, para1);
        [T4, NI4, G4, x4] = IDFPM(fun, 'MAlg1',         3,  8, x0, para1);
        [T5, NI5, G5, x5] = IDFPM(fun, 'TQYAlg3p1',      3, 10, x0, para1);
        [T6, NI6, G6, x6] = IDFPM(fun, 'TVCAlg1',   3, 12, x0, para1);
        [T7, NI7, G7, x7] = IDFPM(fun, 'TAlg3p1',      3, 14, x0, para1);

        % Store NI / T / G (21 columns)
        progress_r = [progress_r; ...
            NI1, T1, G1, ...
            NI2, T2, G2, ...
            NI3, T3, G3, ...
            NI4, T4, G4, ...
            NI5, T5, G5, ...
            NI6, T6, G6, ...
            NI7, T7, G7];

        % Accuracy for all 7 algorithms (using sign for ±1 labels)
        pred1 = sign(A*x1); acc1 = mean(pred1 == b);
        pred2 = sign(A*x2); acc2 = mean(pred2 == b);
        pred3 = sign(A*x3); acc3 = mean(pred3 == b);
        pred4 = sign(A*x4); acc4 = mean(pred4 == b);
        pred5 = sign(A*x5); acc5 = mean(pred5 == b);
        pred6 = sign(A*x6); acc6 = mean(pred6 == b);
        pred7 = sign(A*x7); acc7 = mean(pred7 == b);

        acc_r  = [acc_r;  acc1, acc2, acc3, acc4, acc5, acc6, acc7];
        time_r = [time_r; T1,   T2,   T3,   T4,   T5,   T6,   T7  ];

    end

    % Check and extract best metrics
    if size(progress_r,2) == 21
        TM = min(progress_r, [], 1);  % 1x21 (best over repeats, column-wise)
        % Optionally print them to console
        fprintf(['TM: ', repmat('%g ',1,21), '\n'], TM);
    else
        fprintf('Error: progress_r has incorrect size (%d). Skipping this dataset.\n', size(progress_r,2));
        continue;
    end

    % Per-dataset best accuracy (per algorithm) and fastest training time
    acc_mean          = max(acc_r, [], 1);   % 1x7
    training_time_min = min(time_r, [], 1);  % 1x7

    accuracy_results      = [accuracy_results;      acc_mean];
    training_time_results = [training_time_results; training_time_min];

    % Write one LaTeX-friendly line to file:
    % dataset & m & n & (NI/T/G)x7 & acc1 ... acc7
    fprintf(fid_tex, '%s & %d & %d', dataset{i}, m, n);

    idxNI = [1,4,7,10,13,16,19];
    idxT  = [2,5,8,11,14,17,20];
    idxG  = [3,6,9,12,15,18,21];

    for a = 1:7
        fprintf(fid_tex, ' & %.0f/%.4f/%.2e', TM(idxNI(a)), TM(idxT(a)), TM(idxG(a)));
    end
    for a = 1:7
        fprintf(fid_tex, ' & %.4f', acc_mean(a));
    end
    fprintf(fid_tex, ' \\\\ \r\n');
end

% Close file
fclose(fid_tex);

% Plot training time (grouped bars: datasets x 7 algs)
figure;
bar(training_time_results, 'grouped');
title('Training Time Comparison');
xlabel('Dataset');
ylabel('Training Time (s)');
set(gca,'XTick',1:length(dataset),'XTickLabel',dataset,'XTickLabelRotation',45);
legend({'PASY Alg. ','CTC Alg.','DCZR Alg.',' M Alg. 3.1','TQY Alg.','TVC Alg.','T Alg.'}, ...
       'Location','best');

% Display accuracy results
disp('Accuracy Results (rows = datasets, cols = 7 algorithms):');
disp(accuracy_results);

% Display training time results
disp('Training Time Results (rows = datasets, cols = 7 algorithms):');
disp(training_time_results);
