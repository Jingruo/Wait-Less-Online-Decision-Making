clear; clc; close all;

m = 1;
earlystop = 0;

num_n = 5;
ntest = 100;

type_list = [1, 2];
for type_num = 1 : length(type_list)
    type = type_list(type_num);
    reg_subgradient = zeros(num_n, ntest);
    reg_two_path_know_mu = zeros(num_n, ntest);
    reg_freq = zeros(num_n, ntest);
    reg_freq2 = zeros(num_n, ntest);
    
    time = zeros(num_n, 4);

    % Initialize runtime storage
    runtime_subgradient = zeros(num_n, ntest);
    runtime_two_path_know_mu = zeros(num_n, ntest); 
    runtime_freq = zeros(num_n, ntest);
    runtime_freq2 = zeros(num_n, ntest); 

    vec_A = cell(num_n, ntest);
    vec_b = cell(num_n, ntest);
    vec_c = cell(num_n, ntest);

    idx = 0;
    for n = ceil(logspace(2, 6, num_n))
        idx = idx + 1;
        for itest = 1 : ntest
        % parfor itest = 1 : ntest  % use for parallel computing
            [A, b, c] = olpgetdata(m, n, type);
            vec_A{idx, itest} = sparse(A);
            vec_b{idx, itest} = b;
            vec_c{idx, itest} = c;
            
            % solve for offline
            [xopt, yopt] = olpgurobi(A, b, c);
            objval = c' * xopt;

            % Strategy 1: subgradient
            tstart = tic;
            [x1, y1_list] = olpsubgrad(A, b, c, earlystop);

            % run time needs to be tracked without parallel computation
            runtime_subgradient(idx, itest) = toc(tstart);
            time(idx, 1) = time(idx, 1) + toc(tstart);

            reg_subgradient(idx, itest) = norm(max(A * x1 - b, 0)) + objval - c' * x1;
            
            % Strategy 2: two path first-order
            mu = 1;
            tstart = tic;
            [x2, y2_list] = olptwopath_grad(A, b, c, mu, earlystop);
            runtime_two_path_know_mu(idx, itest) = toc(tstart);
            time(idx, 2) = time(idx, 2) + toc(tstart);
            reg_two_path_know_mu(idx, itest) = norm(max(A * x2 - b, 0)) + objval - c' * x2;
            
            % Strategy 3: main algorithm 1
            mu = 1;
            tstart = tic;
            [x3, y3_list] = olptwopath_freq(A, b, c, mu, earlystop);
            runtime_freq(idx, itest) = toc(tstart);
            time(idx, 3) = time(idx, 3) + toc(tstart);
            reg_freq(idx, itest) = norm(max(A * x3 - b, 0)) + objval - c' * x3;

            % Strategy 4: main algorithm 2
            mu = 1;
            tstart = tic;
            [x4, y4_list] = olptwopath_freq2(A, b, c, mu, earlystop);
            runtime_freq2(idx, itest) = toc(tstart);
            time(idx, 4) = time(idx, 4) + toc(tstart);
            reg_freq2(idx, itest) = norm(max(A * x4 - b, 0)) + objval - c' * x4;

            % Other strategies like pure LP can be adjusted from code in main algorithm 1

        end % End for
    end
end


% plot for log-scale

n_list = logspace(2, 6, num_n);

reg_1 = reg_subgradient;

reg_2 = reg_two_path_know_mu;

reg_3 = reg_freq;

reg_4 = reg_freq2;

reg_1_plot = mean(reg_1, 2); 
reg_2_plot = mean(reg_2, 2);
reg_3_plot = mean(reg_3, 2);
reg_4_plot = mean(reg_4, 2);

loglog(n_list, reg_1_plot, '-', 'LineWidth', 2, 'MarkerSize', 12);

hold on

loglog(n_list, reg_2_plot, '-', 'LineWidth', 2, 'MarkerSize', 12);

loglog(n_list, reg_3_plot, '-', 'LineWidth', 2, 'MarkerSize', 12); 

loglog(n_list, reg_4_plot, '-', 'LineWidth', 2, 'MarkerSize', 12);

loglog(n_list, reg_5_plot, '-', 'LineWidth', 2, 'MarkerSize', 12);
 
legend('offline', 'first-order', 'main 1', 'main 2');


