function [L_hat, P_hat, S_hat, t_hat, ...
    P_track_full, T_calc]= NORST_random(M, T, r, ev_thresh, alpha, K,R,hop)
%Algorithm that implments the NORST random function for subspace tracking
%with missing data or dynamic matrix completion
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
%thresh = ev_thresh / 2;
%[~, r_init] = size(P_init);


[n, t_max] = size(M);

P_init = orth(randn(n, r));

P_hat = P_init;

%T_hat = zeros(n, t_max);
S_hat = zeros(n, t_max);
S_hat_buffer = zeros(n, t_max);
L_hat = zeros(n, t_max);
L_hat_buffer = zeros(n, t_max);
t_hat = [];

k = 0;
cnt = 1;
P_track_full{cnt} = P_hat;
T_calc(cnt) = 1;
ph = 1;     %ph - 0 => detect, 1 => pca
nextbatch_start = 0;
alpha0 = alpha;
permit = 0;
%% Main Algorithm Loop
for ii = 1 : t_max
    %% Estimate support
    Atf.times = @(x) x - (P_hat * (P_hat' * x));
    Atf.trans = @(y) y - (P_hat * (P_hat' * y));
    phi.times = @(x) x - (P_hat * (P_hat' * x));
    y_t = Atf.times(M(:, ii));
%     opts.tol   = 1e-4;
%     opts.print = 0;
%     
%     opts.delta = omega * 2 / 15;
%     
%     x_t_hat_cs = yall1(Atf, y_t, opts);
%     omega = sqrt(M(:, ii)' * M(:, ii) / n);
%     
%     t_hat_temp = find(abs(x_t_hat_cs) > omega);
    %T_hat(t_hat_temp, ii) = 255;
    
    
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
    if (ii == 1)
        x0 = zeros(n,1);
    else
%         x0 = S_hat(:,ii-1);
    end
    
    T_union = find(T(:, ii) == 0);
    tol = 1e-16;
    S_hat(T_union, ii) = ccgls(@Phifun, y_t, T_union, P_hat, ...
        0, tol, 20);
    L_hat(:, ii) = M(:, ii) - S_hat(:, ii);    
    
    %% Subspace update
    
%     if(~mod(ii + 1 , alpha))    
    if(ii == alpha0)
        permit = 1;
    end
    if (ii == nextbatch_start + alpha)
        permit = 1;
    end
    if(permit == 1)                
        idx = nextbatch_start + 1 : nextbatch_start + alpha;
        nextbatch_start = idx(hop);        

        if(idx(end) > t_max)
            idx = nextbatch_start : t_max;
        end
        
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
        else        %%update phase

            P_hat = simpleEVD(L_hat(:, idx), r);
            if( ii==idx(end) && R > 0)
%                     fprintf("fine tuning the subspace at %d\n",ii);
                    idx_buffer = idx;
%                     ctr = [ctr, ii];
%                     idx_buffer = ctr(end - 1) + 1 : ctr(end);
%                     idx_buffer = t_hat(end) : ii;
                for reuse = 1:R
                    for kk = idx_buffer
                        Atf.times = @(x) x - (P_hat * (P_hat' * x));
                        Atf.trans = @(y) y - (P_hat * (P_hat' * y));
                        y_t_buffer = Atf.times(M(:, kk));
                        
%                         x_t_hat_cs_buffer = yall1(Atf, y_t_buffer, opts);                        

%                         x_cs_hat_buffer(:,kk) = x_t_hat_cs_buffer;
                        
%                         omega = 1 * sqrt(M(:, kk)' * M(:, kk) / n);
%                         t_hat_temp_buffer = find(abs(x_t_hat_cs_buffer) > omega);
%                         T_hat_buffer(t_hat_temp_buffer, kk) = 255;
        
%                         T_union_buffer = t_hat_temp_buffer;
                        T_union_buffer = find(T(:,kk) == 0);
                        S_hat_buffer(T_union_buffer, kk) = ccgls(@Phifun, y_t_buffer, ...
                            T_union_buffer, P_hat, ...
                            0, tol, 20);
                        L_hat_buffer(:, kk) = M(:, kk) - S_hat_buffer(:, kk);                        
                    end
%                     permuted_idx = randperm(buffer) + nextbatch_start - alpha;
                    P_hat = simpleEVD(L_hat_buffer(:, idx_buffer), r);
                end
            end
            

%             phi_t = speye(n) - P_hat * P_hat';
            
            k = k + 1;
            if k == K
                ph = 0;                
            end

        end
    
    %% Return subspace
%    if((ii == 1) || ~(mod(ii + 1, alpha)))
        cnt = cnt + 1;
        P_track_full{cnt} = P_hat;        
        T_calc(cnt) = ii;        
        permit = 0;
        
    end
end
end
