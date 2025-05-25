function [x, y_list] = olpsubgrad(A, b, c, earlystop, varargin)
% Use raw subgradient

[m, n] = size(A);

if nargin > 4
    alpha = 1 / n^(varargin{1});
else
    alpha = 1 / n^(0.5);  
end
    

y_list = [];
d = b / n;
y = zeros(m, 1);
x = zeros(n, 1);
bleft = b;

for i = 1:n
    
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
    y = y - alpha * (d - Ai * xi);
    y = max(y, 0);
    y_list = [y_list, y];
    
end % End for

end % End function