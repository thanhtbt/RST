function [U_est,PER] = NORST_ST(M,U_tr,OPT)

%This folder contains the code accompanying pre-print.
%[1] "Nearly Optimal Robust Subspace Tracking", Praneeth Narayanamurthy and Namrata Vaswani, arXiv:1712.06061, 2017.

%%%                          Inputs                         %%%
%%%     M - measurement matrix                              %%%
%%%     ev_thres - threshold for subspace change detection  %%%
%%%     P_init - an initial estimate of the subspace        %%%
%%%     t_train - the dimension of the training data        %%%


%%%                       Algorithm parameters              %%%
%%%     alpha - frame length                                %%%
%%%     K - number of projection PCA steps                  %%%
%%%     omega - threshold for non-zero value in S           %%%
%%% 	K_CS - number of CoSaMP iterations 		    %%%
%%% 	outc - an upper bound on estimate of fraction of outliers per column

%%%                          Outputs                        %%%
%%%     L_hat - estimate of the low rank matrix             %%%
%%%     P_hat - estimate of the subspace in which L lies    %%%
%%%     S_hat - estimate of the sparse signal               %%%
%%%     t_hat - estimate of subspace change times           %%%

%% Initializations

[n, t_max] = size(M);
[r] = size(U_tr,2);

if isfield(OPT,'K'),
    K = OPT.K;
else
    K = 5;
end
if isfield(OPT,'omega'),
    omega = OPT.omega;
else
    omega = 5;
end
if isfield(OPT,'ev_thresh'),
   ev_thresh  = OPT.ev_thresh;
else
    ev_thresh = 2e-3;
end

if isfield(OPT,'t_train'),
   t_train  = OPT.t_train;
else
    t_train = 200;
end
if isfield(OPT,'outc'),
   outc  = OPT.outc;
else
    outc = ceil(0.05 * n); 
end
if isfield(OPT,'alpha'),
    alpha  = OPT.alpha
else
    alpha = 200;
end
if isfield(OPT,'K_CS'),
   K_CS = OPT.K_CS; 
else
   K_CS = 5;
end

M_train = M(:,1:t_train);
P_init  = orth(M_train);
P_init  = P_init(:,1:r);

PER.EV  = zeros(1,t_max);
PER.SEP = zeros(1,t_max);
PER.SE  = zeros(1,t_max);
PER.Angle = zeros(1,t_max);


[~, r_init] = size(P_init);
P_hat_old   = P_init;
P_hat_new   = [];
P_hat       = [P_hat_old, P_hat_new];


T_hat = zeros(n, t_max);
S_hat = zeros(n, t_max);
L_hat = zeros(n, t_max);
t_hat = [];

% L_hat(:, 1 : t_train) = M(:, 1 : t_train);
k = 0;
cnt = 1;
phi_t = (eye(n) - P_hat * P_hat');
% In = eye(n);
ph = 1;     %ph - 0 => detect, 1 => ppca
opts.delta=0.4;

%% Main Algorithm Loop
for ii = 1 : t_max
    %% Estimate support
    Atf.times = @(x) x - (P_hat * (P_hat' * x));
    Atf.trans = @(y) y - (P_hat * (P_hat' * y));
    phi.times = @(x) x - (P_hat_old * (P_hat_old' * x));
    y_t = Atf.times(M(:, ii));
    opts.tol   = 1e-4;
    opts.print = 0;
    
    opts.delta = omega * 2 / 15;
    
    x_t_hat_cs = yall1(Atf, y_t, opts);
    %omega = sqrt(M(:, ii)' * M(:, ii) / n);
    
    t_hat_temp = find(abs(x_t_hat_cs) > omega);
    T_hat(t_hat_temp, ii) = 255;
    
    %LS.times = @(x) phi(:, t_hat_temp) * x;
    %LS.trans = @(y) phi(:, t_hat_temp)' * x;
    
    %y_t = M(:, ii) - (P_hat * (P_hat' * M(:, ii)));
    %DecayRate = 0.9; %other values work, may make it slower
    %x_t_hat_cs = cosamp_cgls(phi_t, ...
    %    y_t, outc, DecayRate, K_CS, 1e-6);
    %t_hat_temp = find(abs(x_t_hat_cs) > omega);
    %     T_hat(t_hat_temp, ii) = 1;
    
    %% Estimate signal components
    % %         [S_hat(t_hat_temp, ii), ~] = ...
    % %             lsqr(phi_t(:, t_hat_temp), y_t, 1e-6, 50);
    %     S_hat(t_hat_temp, ii) = phi_t(:, t_hat_temp) \ y_t;
    S_hat(t_hat_temp, ii) = cgls(phi_t(:, t_hat_temp), y_t, ...
        0, 1e-10, 10);
    L_hat(:, ii) = M(:, ii) - S_hat(:, ii);
    
    
    %% Subspace update
    if(~mod(ii + 1 , alpha))
        u = (ii + 1) / alpha;
        idx = (u-1) * alpha + 1 : u * alpha ;
        L_temp = L_hat(:, idx);
        
        %MM = L_temp - (P_hat_old *(P_hat_old' * L_temp));
        MM = phi.times(L_temp);
        
        if(~ph)     %%detect phase
            %             phi_t = eye(n) - P_hat * P_hat';
            if(svds(MM, 1) >= sqrt(alpha * ev_thresh))
                ph = 1;
                t_hat = [t_hat, ii];
                k = 0;
            end
        else        %%ppca phase
            P_hat = simpleEVD((L_hat(:, max(1, ii - alpha + 1) : ii)), r_init);
            phi_t = speye(n) - P_hat * P_hat';
            k = k + 1;
            
            if(k==K + 1)
                P_hat_old = P_hat;
                k = 1;
                ph = 0;
                phi_t = speye(n) - P_hat * P_hat';
            end
        end
    %end
    
    
%     %% Return subspace
% %    if((ii == 1) || ~(mod(ii + 1, alpha)))
%         P_track_full{cnt} = P_hat;
%         %P_track_new{cnt} = P_hat_new;
%         T_calc(cnt) = ii;
%         cnt = cnt + 1;
    end
    
    PER.SEP(1,ii) = sub_est_per(U_tr , P_hat, 'SEP');
    PER.Angle(1,ii)  = sub_est_per(U_tr , P_hat, 'Angle');
    PER.SE(1,ii) = sub_est_per(U_tr , P_hat, 'SE');
    PER.EV(1,ii) = sub_est_per(U_tr , P_hat, 'EV');
end
U_est = P_hat;
end
