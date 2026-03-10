clc;
clear all;
close all;

% Set random seed
seed = 97006855;
ss = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(ss);

% Set parameters
ITR_max = 2000;
lambda = 1e-1;
para1.Itr_max = ITR_max;
para1.gamma = 1;
para1.sigma = 1;
para1.tau = 0.5;
para1.rho = 1.5;
para2.Itr_max = ITR_max;
para2.gamma = 1;
para2.sigma = 0.001;
para2.tau = 0.5;
para2.rho = 1.8;

% Dataset
dataset = {'a1a.t', 'a2a.t', 'a3a.t', 'a4a.t', 'a5a.t', 'a6a.t', 'a7a.t', 'a8a.t', 'a9a.t', 'colon-cancer'};
fid_tex = fopen('Rel.txt', 'w');
accuracy_results = [];
time_results = [];

for i = 1:length(dataset)
    [b, A] = libsvmread(dataset{i});
    [m, n] = size(A);
    fprintf('name=%s, m=%d, n=%d, lambda=%.2f\n', dataset{i}, m, n, lambda);
    fun = @(x) lr_loss(A, b, m, x, lambda);
    progress_r = [];
    acc_r = [];
    time_r = [];
    for repeat = 1:5
        % Set the initial point
        x0 = 4 * (rand(n, 1) - 0.5);
        % Start comparison
        [T1, NI1, G1, x1] = IDFPM(fun, 'DIALG', 3, 2, x0, para2);
        [T2, NI2, G2, x2] = IDFPM(fun, 'SIALG', 3, 4, x0, para1);
        [T3, NI3, G3, x3] = IDFPM(fun, 'WIALG', 3, 6, x0, para1);
      
        progress_r = [progress_r; NI1, T1, G1, NI2, T2, G2, NI3, T3, G3];
        
        
        % Calculate accuracy
        pred1 = A * x1 >= 0;
        pred2 = A * x2 >= 0;
        pred3 = A * x3 >= 0;
        acc1 = sum(pred1 == b) / length(b);
        acc2 = sum(pred2 == b) / length(b);
        acc3 = sum(pred3 == b) / length(b);
        acc_r = [acc_r; acc1, acc2, acc3];
        time_r = [time_r; T1, T2, T3];
    end
    TM = min(progress_r);
    fprintf('%d, %d, %d, %d, %d, %d, %d, %d, %d, \n', TM(1), TM(2), TM(3), TM(4), TM(5), TM(6), TM(7), TM(8), TM(9));
    acc_mean = mean(acc_r);
    time_mean = mean(time_r);
    accuracy_results = [accuracy_results; acc_mean];
    time_results = [time_results; time_mean];

    fprintf(fid_tex, '%s & %d & %d & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e & %.4f & %.4f & %.4f\\\\ \r\n',...
                dataset{i}, m, n, TM(1), TM(2), TM(3), TM(4), TM(5), TM(6), TM(7), TM(8), TM(9), acc_mean(1), acc_mean(2), acc_mean(3));
end

% Close file
fclose(fid_tex);

% Plot training time
figure;
bar(time_results);
title('Training Time Comparison');
xlabel('Dataset');
ylabel('Training Time (s)');
legend('DIALG', 'SIALG', 'WIALG');

% Display accuracy results
disp('Accuracy Results:');
disp(accuracy_results);
