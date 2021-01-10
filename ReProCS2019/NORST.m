function [L_hat, P_hat, S_hat,T_hat,...
          t_hat,P_track_full, T_calc] = NORST(M,T_obs,P_init,ev_thresh,...
                                              alpha, K, omega, tol)
% This MATLAB function implements the code accompanying pre-print for
% simulated data low-rank matrix recovery:
% [1] "Subspace Tracking from Missing and Outlier Corrupted Data"
%      P. Narayanamurthy, V. Daneshpajooh, N. Vaswani
%      arXiv:1810.03051v1 [cs.LG] 6 Oct 2018
%  
% Input:
% M = measurement matrix
% T_obs = observed entries' support
% P_init = an initial estimate of the subspace
% 
%   Algorithm parameters:
%   ev_thres = threshold for subspace change detection
%   alpha = frame length
%   K = number of projected PCA steps
%   omega = the threshold for non-zero values in S
%   tol = tolerance for cgls function (conjugate gradient LS)
% 
% Output:
% L_hat = estimate of the low rank matrix
% P_hat = estimate of the subspace in which L lies
% S_hat = estimate of the sparse matrix
% t_hat = estimate of subspace change times

%% Initializations
[~, r] = size(P_init);
[n, t_max] = size(M);

P_hat_old = P_init;
P_hat_new = [];
P_hat = [P_hat_old, P_hat_new];

T_hat = zeros(n, t_max);
S_hat = zeros(n, t_max);
L_hat = zeros(n, t_max);

t_hat = [];

k = 0;
cnt = 1;
ph = 1;     %ph : 0 => detect, 1 => ppca

phi_t = (eye(n) - P_hat * P_hat');

opts.delta=0.4;
opts.print = 0;

%% Main Algorithm Loop
for ii = 1 : t_max
    %% Estimate support
    Atf.times = @(x) x - (P_hat * (P_hat' * x));
    Atf.trans = @(y) y - (P_hat * (P_hat' * y));
    phi.times = @(x) x - (P_hat_old * (P_hat_old' * x));
    y_t = Atf.times(M(:, ii));
    opts.tol   = 1e-4;   
    
    weights = ones(n,1); 
    lambda1 = 0;
    weights(find(T_obs(:,ii) == 0)) = lambda1;
    
    opts.weights = weights;
    opts.delta = omega * 2 / 15;
    
    x_t_hat_cs = yall1(Atf, y_t, opts);    
    x_cs_hat(:,ii) = x_t_hat_cs;
    
    t_hat_temp = find(abs(x_t_hat_cs) > omega);
    T_hat(t_hat_temp, ii) = 255;
    
    %% Estimate signal components    
    T_union = unique([t_hat_temp;find(T_obs(:,ii) == 0)]);
    
    S_hat(T_union, ii) = cgls(phi_t(:, T_union), y_t, 0, tol, 20);
    L_hat(:, ii) = M(:, ii) - S_hat(:, ii);    
    
    %% Subspace update
    if(~mod(ii + 1 , alpha))
        u = (ii + 1) / alpha;
        idx = (u-1) * alpha + 1 : u * alpha ;
        
        L_temp = L_hat(:, idx);       
        MM = phi.times(L_temp);        
        if(~ph)     %%detect phase            
            if(svds(MM, 1) >= sqrt(alpha * ev_thresh))
                ph = 1;
                t_hat = [t_hat, ii];
                k = 0;
            end
        else        %%pca phase
            P_hat = simpleEVD((L_hat(:, max(1, ii - alpha + 1) : ii)), r);
            phi_t = speye(n) - P_hat * P_hat';
            k = k + 1;            
            if(k == K + 1)
                P_hat_old = P_hat;
                k = 1;
                ph = 0;
                phi_t = speye(n) - P_hat * P_hat';
            end
        end    
        %% Return subspace
        P_track_full{cnt} = P_hat;        
        T_calc(cnt) = ii;
        cnt = cnt + 1;
    end
end
end