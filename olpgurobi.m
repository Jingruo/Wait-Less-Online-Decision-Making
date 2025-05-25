function [x, y] = olpgurobi(A, b, c)

[~, n] = size(A);

model.A = sparse(A);
model.rhs = b;
model.obj = c;
model.sense = '<';
model.modelsense = 'max';
model.lb = zeros(n, 1);
model.ub = ones(n, 1);
param.Presolve = 0;
param.Method = 1;
param.OutputFlag = 0;

sol = gurobi(model, param);
x = sol.x;
y = sol.pi;

end % End function