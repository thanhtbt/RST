% Copyright (c) 2017 Mitsubishi Electric Research Laboratories (MERL). All rights reserved.
% 
% The software, documentation and/or data in this file is provided on an "as is" 
% basis, and MERL has no obligations to provide maintenance, support, 
% updates, enhancements or modifications. MERL specifically disclaims any 
% warranties, including, but not limited to, the implied warranties of 
% merchantability and fitness for any particular purpose. In no event shall 
% MERL be liable to any party for direct, indirect, special, incidental, or 
% consequential damages, including lost profits, arising out of the use of 
% this software and its documentation, even if MERL has been advised of the 
% possibility of such damages.
% 
% As more fully described in the license agreement that was required in order
% to download this software, documentation and/or data, permission to use, 
% copy and modify this software without fee is granted, but only for 
% educational, research and non-commercial purposes.

function [L, R, S, ParamsNew] = ROSETA(A, b, L_prev, r, Params, lambda)
%
% Solves the subspace tracking problem using the augmented Lagrangian method:
%
%   min_{L, R, S}:   lambda||S||_1 + ||b - A(L*R' + S)||_2^2 +  \|L_prev - L\|_H^2
%   
%
% Input:
% A - mask of current frame (can be empty);
% b - current frame;
% L_prev - previous estimate of the subspace
% r - rank of L;
% Params: a struct saving parameters
%   Params.const - sets how fast the stepsize can change;
%   Params.temp - the parameter to compute the stepsize mu =
%   Params.const/(1+Params.temp). This will not be used to determine the current
%   stepsize, but will be used in the computation of the new stepsize. This
%   is kept in the Params mainly for debugging purposes.
%   Params.HIGH_TEMP - maximum value of hyper parameter temp
%   Params.mu - current stepsize;
%   Params.mufloor - minimum value of stepsize mu
%   Params.gd - the gradient of last frame
%
% Output:
% L - estimated subspace;
% R - subspace coefficients;
% S - sparse outliers of current frame;
% ParamsNew - updated parameters;

%
% Written by Hassan Mansour (mansour@merl.com), Xin Jiang (chlorisjiang@gmail.com)
% Copyright MERL 2017

mu = Params.mu;
val = Params.val;
ParamsNew = Params;
ParamsNew.norminner = 0;


HIGH_TEMP = Params.HIGH_TEMP;
ParamsNew.HIGH_TEMP = HIGH_TEMP;



mufloor = Params.mufloor;
ParamsNew.mufloor = mufloor;

m = length(b);

if ~exist('mu', 'var')
    mu = 1.2/norm(b(:));
end


B = b;

if isempty(A)
    mask = ones(m,1);
else
    mask = A>0;
end

% initialize the variables
S = zeros(m,1);
E = zeros(m,1);

if size(L_prev, 2) == r %rank(L_prev) == r
    U = L_prev;
    
    %% Find coefficient and sparse updates
    
    options.tol                 = 1e-6;
    options.maxIter             = 100;
    options.mask                = mask;
    
    [  S(mask), alpha] = SubspaceFitting( U(mask,:), B(mask), options, lambda ); 
    E = -U*alpha;
    E(mask) = 0;
    
    
    %% find low-rank updates
    
    
    G = (B - S - E - L_prev*alpha) * alpha';
    Dec = (G - (G*alpha)*alpha'/(1 + alpha'*alpha));
    
    % update subspace
    U = L_prev + Dec/mu;
    
    %% deciding the next stepsize mu
    % compute the current gradient
    if size(U,2) == r
        gd = Dec;
        ParamsNew.gd = gd;
        ParamsNew.inner = 0;
        ParamsNew.norminner = norm(gd,'fro');
    end
    
    if isfield(Params,'gd') % if previous gradient is available
        % compute the trace of the inner product
        ParamsNew.inner = trace(Params.gd' * ParamsNew.gd);
        
        % normalize the inner product so it is not too large/small
        norm_inner = norm(Params.gd,'fro')*norm(ParamsNew.gd,'fro');
        ParamsNew.norminner = norm_inner;
        % two tuning parameters
        const = Params.const;
        val = Params.val; % decides how fast the stepsize (\mu) can change
        
        Params.val = val;

        ParamsNew.inner = ParamsNew.inner / norm_inner;

        % compute the hyper-parameter that decides the value of \mu
        ParamsNew.temp = min(max(Params.temp + sigmoid(ParamsNew.inner,val), const*1e-1), HIGH_TEMP);

        % compute \mu for next frame
        mu_new = (const)/(1+ParamsNew.temp);
        
    end
    
    % if no update is computed (previous gradient information unavailable), keep \mu constant
    if ~exist('mu_new','var')
        mu_new = mu;
    end
    % update the new \mu
    ParamsNew.mu = mu_new;
    ParamsNew.val = val;
    
    X = U*alpha;
   
    %% Print progress
    res = norm(b(:) - A.*(X(:) + S(:) + E(:)));
    
    
else % if rank(U_prev) != r, then enter the training stage
    
    if isempty(L_prev) % no initialization
        alpha = norm(B);
        U = B./alpha;
    else
        L_prev = orth(L_prev);
        u = (B - L_prev * (L_prev' * B))./norm(B - L_prev * (L_prev' * B));
        U = cat(2, L_prev, u);
        alpha = pinv(U'*U)*U'*B;
        
    end
end


S = A.*S(:);
L = U;
R = alpha';


end

function fval = sigmoid(x, val)
FMIN = -val; FMAX = val;
omega = 0.1;

fval = FMIN + (FMAX - FMIN)/(1 - (FMAX/FMIN)*exp(-x/omega));
end

function [  s, a ] = SubspaceFitting( U, b, options, lambda )
%
% [ s, a] = SubspaceFitting( U, b, options, mu, lambda )
% 
% Solve the subspace fitting with sparse outliers problem using Lasso
% 
%   minimize_{a, s}     1/2||Ua + s - b||_2^2 + \lambda||s||_1 
%
%   
% Written by Xin Jiang (chlorisjiang@gmail.com), Hassan Mansour (mansour@merl.com), 
% Copyright MERL, 2017

maxIter = options.maxIter;
tol = options.tol;

% initialize variables
[m, n] = size(U);

a = zeros(n,1);
s = zeros(m,1);



%% Solver
pinvU = pinv(U);
for iter = 1:maxIter
    
    % update subspace coefficients
    a = pinvU * (b -s);
    Ua = U*a;

    % restrict magnitude to observed measurements
    Ua(Ua>max(b(:))) = max(b(:));
    Ua(Ua<min(b(:))) = min(b(:));

    % update outlier vector
    s = b-Ua;
    s = sign(s).*max(0, abs(s) - lambda);

    % update
    res = Ua + s - b;
    
    if (norm(res) < tol)
        break;
    end
   
end

end

