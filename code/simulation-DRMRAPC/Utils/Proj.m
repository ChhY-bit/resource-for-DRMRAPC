function res = Proj(x,y,f,grad_f)
if nargin == 3
    n = length(x); %#ok<NASGU>
    syms x_sym [n,1]
    grad = jacobian(f(x_sym),x_sym);
    grad_f = matlabFunction(grad, 'Vars', {x_sym});
end
    g = grad_f(x)';
if ~iscolumn(g)
    g = g';
end
if ~iscolumn(y)
    y = y';
end
    if y'*g > 0 && f(x) > 0
        res = y - g*g'/norm(g)^2*y*f(x);
    else
        res = y;
    end
end