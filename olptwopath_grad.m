function [x, y_list, tau] = olptwopath_grad(A, b, c, mu, earlystop)

% Use two paths, know mu

[m, n] = size(A);
y_list = [];

d = b / n;
n1 = ceil(n^(2/3));
y = zeros(m, 1);
y_second = y;
x = zeros(n, 1);
bleft = b;
alpha = 1 / (n)^(1/3);
tau = 0;

for i = 1:n1
    
    ci = c(i);
    Ai = A(:, i);
    
    % first path: subgradient, and make decision
    if ci > Ai' * y
        bleft = bleft - Ai;
        if min(bleft) < 0 && earlystop
            break;
        end % End if
        xi = 1;
    else
        xi = 0;
    end % End if
    
    x(i) = xi;
    y = y - alpha * (d - Ai * xi);
    y = max(y, 0);
    y_list = [y_list, y];
    
    % second path: know mu, only learn, no decision
    alpha_second = 2 / (mu * (i+1));
    
    if ci > Ai' * y_second
        xi_null = 1;
    else
        xi_null = 0;
    end % End if
    
    y_second = y_second - alpha_second * (d - Ai * xi_null);
    y_second = max(y_second, 0);
   
end % End for

alpha = 1 / (n)^(2/3);

y = y_second;

for i = (n1+1):n
    
    ci = c(i);
    Ai = A(:, i);
    
    if ci > Ai' * y
        bleft = bleft - Ai;
        if min(bleft) < 0 && earlystop
            break;
        end % End if
        xi = 1;
    else
        xi = 0;
    end % End if

    x(i) = xi;
    
    x(i) = xi;
    y = y - alpha * (d - Ai * xi);
    y = max(y, 0);
    
    y_list = [y_list, y];
    
end % End for

end % End function