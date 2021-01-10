function P_hat = simpleEVD(X, r)
%   This MATLAB code implements the Simple-EVD algorithm and returns a
%   basis matrix for the new subspace
%
%   Input:
%   X = emperical covariance data matrix (m x m)
%   r = target rank of output
%
%   Output:
%   P_hat = basis matrix for output (m x r)

[P_hat, ~] = svds(X, r);
end
