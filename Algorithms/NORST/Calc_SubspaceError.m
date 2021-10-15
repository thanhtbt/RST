function[SE_err] = Calc_SubspaceError(A, B)
%   This MATLAB code calculates the subspace error of two subspaces
%
%   Input:
%   A,B = matrices of size (m x n)
%
%   Output:
%   SE_err = subspace error, i.e, sin(\theta) where \theta is the angle b/w
%            the subspaces where A and B lie in.

if(~isempty(A) && ~isempty(B))
    [QA, ~] = qr(A);
    [QB, ~] = qr(B);
    
    ra = rank(A);
    rb = rank(B);
    
    [m,~] = size(A);
    SE_err = norm((eye(m) - QA(:, 1 : ra)  * QA(:, 1 : ra)') ...
        * QB(:, 1 : rb), 2);
else
    error('Error: empty input argument!\n');
end

end
