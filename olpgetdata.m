function [A, b, c] = olpgetdatarebuttal(m, n, type)

% Produces random data for Online Linear Programming (OLP)
%                  experiments in resource allocation.
%
%   [A, b, c] = olpgetdata(m, n, dataType)
%   returns:
%       - A: an m-by-n matrix of resource usage coefficients,
%       - b: an m-by-1 resource capacity vector,
%       - c: an n-by-1 revenue/benefit vector,
%   given:
%       - m: number of resource constraints,
%       - n: number of items/customers,
%       - dataType: an integer (1 through 10) determining the distribution
%         from which A, b, and c are sampled.
%
%   Example usage:
%       [A, b, c] = olpgetdata(5, 100, 3);

if type == 1
    A = rand(m, n) * 2;
    b = n * (rand(m, 1) + 1) / 3;
    c = rand(n, 1) * 2;
elseif type == 2
    A = randn(m, n) + 1;
    b = n * (rand(m, 1) + 1) / 3;
    c = sum(A, 1)' - rand(n, 1) * m;
elseif type == 3
    A = 1 + trnd(1, m, n);
    A = min(A, 10);
    A = max(A, -10);
    b = n * (rand(m, 1) + 1) / 3;
    c = sum(A, 1)' - rand(n, 1) * m;    
elseif type == 4
    A = ones(1, n);
    c = rand(n, 1);
    b = n * (rand(1, 1) + 1) / 3;
elseif type == 5
    A = randn(m, n) + 0.5;
    c = sum(A, 1)';
    dlist = [0.2; 0.3];
    b = n * dlist(randi(1:2, m, 1));
elseif type == 6
    A = abs(rand(m, n) * 1.5 - 0.5);
    c = rand(n, 1) * 10;
    b = n * ones(m, 1) * 0.25;
elseif type == 7
    A = ones(1, n);
    c = rand(n, 1);
    b = n / 2;
elseif type == 8
    A = rand(m, n);
    c = rand(1, n);
    selected_indices = randperm(n, 10);
    selected_column = A(:, selected_indices);
    selected_revenue = c(selected_indices);
    expanded_column = repmat(selected_column, 1, n);
    expanded_revenue = repmat(selected_revenue, 1, n);
    idx_set = randperm(n);
    A = expanded_column(:, idx_set);
    c = expanded_revenue(idx_set)';
    b = n * (rand(m, 1) + 1) / 3;
elseif type == 9
    A = 1 + rand(m, n) * 5;
    b = n * (rand(m, 1) + 1) / 3;
    c = rand(n, 1) * 3;
elseif type == 10
    A = betarnd(1, 8, m, n) + 5;
    b = n * (rand(m, 1) + 1) / 3;
    c = rand(n, 1) * 5;
    

end % End if 

end % End function