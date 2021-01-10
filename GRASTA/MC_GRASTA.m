function [X_re,SEP] = MC_GRASTA(X,r,W_ext)
%  GRASTA (Grassmannian Robust Adaptive Subspace Tracking Algorithm) robust 
%  matrix completion code   
%  by Jun He and Laura Balzano, Sept. 2011.
%
%   Online Robust Subspace Tracking from Partial Information
%       http://arxiv.org/abs/1109.3827
%
% Inputs: 
%
%       (I,J,S) index the known entries across the entire data set X. So we
%       know that for all k, the true value of X(I(k),J(k)) = S(k)
%
%       numr = number of rows
%       numc = number of columns
%           NOTE: you should make sure that numr<numc.  Otherwise, use the
%           transpose of X
%
% Outputs: 
%       U and R such that UR' approximates X.
%
%
% Form some sparse matrices for easier matlab indexing

[numr, numc] = size(X);

%=======================GRASTA parameters==================================
OPTIONS.QUIET               = 1;     % suppress the debug information

OPTIONS.MAX_LEVEL           = 20;    % For multi-level step-size,
OPTIONS.MAX_MU              = 15;    % For multi-level step-size
OPTIONS.MIN_MU              = 1;     % For multi-level step-size

OPTIONS.DIM_M               = numr;  % your data's ambient dimension
OPTIONS.RANK                = r; % give your estimated rank

OPTIONS.ITER_MIN            = 200;    % the min iteration allowed for ADMM at the beginning
OPTIONS.ITER_MAX            = 200;    % the max iteration allowed for ADMM
OPTIONS.rho                 = 2;   % ADMM penalty parameter for acclerated convergence
OPTIONS.TOL                 = 1e-8;   % ADMM convergence tolerance

OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.
                                     
CONVERGE_LEVLE              = 20;    % If status.level >= CONVERGE_LEVLE, robust mc converges

OPTS  = struct();   % initial a empty struct for OPTS
D = zeros(1);       % U_hat will be initialized in GRASTA
status = struct();  % maintain GRASTA running status
status.init  = 0;   % will be set 1 once GRASTA start working
%==========================================================================


% Subspace Update
for k = 1:numc,
    x = X(:,k);
    idx = find(x); 
    %====================GRASTA algorithm==================================
    [D, status, OPTS] = grasta_stream(x(idx), idx, D, status, OPTIONS, OPTS);
    if status.level >= CONVERGE_LEVLE,
       break;
    end
end
Wt = D;
SEP = abs(trace(Wt'*(eye(numr)-W_ext*W_ext')*Wt)/trace(Wt'*(W_ext*W_ext')*Wt));
% Reconstruction

OPTS2 = OPTS;
R = zeros(numc,OPTIONS.RANK);
% Outliers = zeros(numc,numr);

for k=1:numc,
    x = X(:,k);
    idx = find(x); 
    x_Omega = x(idx);
    
    if length(idx) < OPTIONS.RANK * 1,
        continue;
    end
        
    D_Omega = D(idx,:);
    
    if OPTIONS.USE_MEX,
        [s, w, ~] = mex_srp(D_Omega, x_Omega, OPTS2);
    else
        [s, w, ~] = admm_srp(D_Omega, x_Omega, OPTS2);
    end

    R(k,:) = w';
   
end
X_re = D*R';
end