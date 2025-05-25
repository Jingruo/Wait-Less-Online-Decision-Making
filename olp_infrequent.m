m = 5;
n = 50;
A = rand(m, n) * 2;
b = (rand(m, 1) + 1) / 3;
r = rand(n, 1) * 2;

sA = sparse(A);
Ts = [1000]';
ls = length(Ts);
lams = ones(n, 1) / n;

argmax_opt_array = zeros(ls, 3);
hindsight_opt_array = zeros(ls, 1);
Alpha = 0.95;
Beta = 0.95;


for l = 1:ls
    T = Ts(l);
    % b = bs(l);
    fprintf('========== Time Horizon %d ========== \n', T)
    simulate_iter = 100;
    hindsight_opt = 0;
    argmax_opt = 0;
    empirical_lams = zeros(n, 1);

    
    time_set = [ceil(T/2)];
    % construct resolving time set 
    for cursor = 1: ceil(log(log(T)/log(3))/log(1/Alpha))
        time_set = [time_set, ceil(T^(Alpha^cursor))];
    end
    time_set = unique(time_set);
    for cursor = 1: ceil(log(log(T)/log(3))/log(1/Beta))
        time_set = [time_set, ceil(T - T^(Beta^cursor))];
    end
    % time_set = [floor(T^(0.51^2)), floor(T/2), floor(T-T^(0.51))];
    
    
    time_set = unique(time_set);
    t_num = length(time_set);
    argmax_time = 0;

    for iter = 1:simulate_iter
        consumption = zeros(m, 1);
        revenue = 0;
        realized_arrival = zeros(n, 1);
        cursor = 1;
        u = zeros(n, 1);
        d = zeros(n, 1);
        
        for t = 1:T
            arrival = mnrnd(1, lams)';
            ind = find(arrival);
            
            %% Argmax Algorithm
            arg = tic;
            if cursor <= t_num
                if t == time_set(cursor)
                    empirical_lams = realized_arrival/max(t-1, 1);
                    [u, ~] = solve_obj(m, n, r, A, b*T-consumption, T-t+1, empirical_lams, 'fluid', 0);
                    model.A = sA;
                    model.obj = r;
                    model.modelsense = 'Max';
                    model.rhs = b*T-consumption;
                    model.lb = zeros(n, 1);
                    model.ub = empirical_lams*(T-t+1);
                    params.outputflag = 0;
                    result = gurobi(model, params);
                    u = result.x;
                    % u = [min(empirical_lams(1)*(T-t+1), b*T-consumption), 0];
                    % u(2) = min(b*T-consumption - u(1), empirical_lams(2)*(T-t+1));
                    d = empirical_lams*(T-t+1);
                    cursor = cursor + 1;
                end
            end

            if consumption + A(:, ind) <= b*T
                if u(ind) >= d(ind)/2
                    % accept
                    revenue = revenue + r(ind);
                    consumption = consumption + A(:, ind);
                    u(ind) = u(ind) - 1;
                end
            end
            d(ind) = d(ind) - 1;
            time_arg = toc(arg);
            argmax_time = argmax_time + time_arg;

            realized_arrival = realized_arrival + arrival;
            
        end

        argmax_opt = argmax_opt + revenue;
        
        %% Hindsight Opt
        % sol = [min(realized_arrival(1), b*T), 0];
        % sol(2) = min(b*T - sol(1), realized_arrival(2));
        % res = sol(1)*r(1) + sol(2)*r(2);
        [~, res] = solve_obj(m, n, r, A, b*T, T, lams, 'hindsight', realized_arrival);
        hindsight_opt = hindsight_opt + res;

        if mod(iter, 10) == 0
            fprintf('Iteration %d/%d \n', iter, simulate_iter)
        end
    end
    
    hindsight_opt_array(l) = hindsight_opt/simulate_iter;
    
    argmax_opt_array(l, 1) = (hindsight_opt-argmax_opt)/simulate_iter;
    argmax_opt_array(l, 2) = argmax_time/simulate_iter;
    argmax_opt_array(l, 3) = t_num;
    
end

data = Ts;
% data = bs';
data(:, 2) = round(argmax_opt_array(:, 1), 1);
data(:, 3) = round(argmax_opt_array(:, 2), 3);
data(:, 4) = round(argmax_opt_array(:, 3), 1);
Table = array2table(data);
Table.Properties.VariableNames(1:4) = {'T', 'argmax_regret', 'argmax_time', 'argmax_LP_num'};


function [sol, res] = solve_obj(m, n, r, A, b, T, lams, type, realized_lam)
% nA = [A; eye(n)];

model.A = sparse(A);
model.obj = r;
model.modelsense = 'Max';
model.rhs = b;
model.lb = zeros(n, 1);
% model.ub = zeros(n, 1);
if strcmp(type, 'fluid')
    model.ub = lams*T;
elseif strcmp(type, 'hindsight')
    model.ub = realized_lam;
end

params.outputflag = 0;

result = gurobi(model, params);

res = result.objval;
sol = result.x;
end