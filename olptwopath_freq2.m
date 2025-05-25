function [x, y_list] = olptwopath_freq2(A, b, c, mu, earlystop)

% Implementation of Main Algorithm 2

[m, n] = size(A);
y_list = [];

d = b / n;
y = zeros(m, 1);
x = zeros(n, 1);
bleft = b;
freq = ceil(n^(1/3));
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
    alpha = 1 / (i + 1)^(2/3);
    y = y - alpha * (d - Ai * xi);  % Subgradient step
    y = max(y, 0);  % Ensure y stays non-negative
    y_list = [y_list, y];

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
    
    else

        alpha = 1 / (n)^(2/3); 
        y = y - alpha * (d - Ai * xi); 
        y = max(y, 0);

    end % End if

    y_list = [y_list, y];

end % End for

end % End function


