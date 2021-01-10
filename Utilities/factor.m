function [L U] = factor(A, rho)
    [m, n] = size(A);
    L = chol(A'*A + rho*speye(n),'lower');
    L = sparse(L);
    U = sparse(L');
end