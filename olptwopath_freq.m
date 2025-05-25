function [x, y_list, tau] = olptwopath_freq(A, b, c, mu, earlystop)

% Implementation of Main Algorithm 1

[m, n] = size(A);
y_list = [];

d = b / n;
y = zeros(m, 1);
x = zeros(n, 1);
bleft = b;
freq = ceil(n^(1/2));
alpha = 1 / freq^(1/2);

for i = 1:n
    ci = c(i);
    Ai = A(:, i);

    % Decision-making via subgradient method
    if ci > Ai' * y
        bleft = bleft - Ai;
        xi = 1;
    else
        xi = 0;
    end

    if min(bleft) <= 0 && earlystop
        break;  % Stop early if resources are depleted
    end

    x(i) = xi;

    % Subgradient update for y
    if i <= freq
        alpha = 2 / (i + 1)^(2/3);
        y = y - alpha * (d - Ai * xi);  % Subgradient step
        y = max(y, 0);
    end

    if i >= n - freq
        alpha = 1 / freq^(2/3);
        y = y - alpha * (d - Ai * xi);  % Subgradient step
        y = max(y, 0);
    end

    % Learning update using LP-based method every 'freq' iterations
    if mod(i, freq) == 0 && min(bleft) > 0
        model.A = sparse(A(:, 1:i));
        model.rhs = bleft*i/(n-i);
        model.obj = c(1:i);
        model.sense = '<';
        model.modelsense = 'max';
        model.lb = zeros(i, 1);
        model.ub = ones(i, 1);
        param.Presolve = 1;
        param.OutputFlag = 0;
    
        sol = gurobi(model, param);
        y = sol.pi;

        y = max(y, 0);

    end % End if

    y_list = [y_list, y];

end % End for

end % End function


