function [U_est,PER] = NORST_ST(X,U_tr,Indicator_full)


[n, t_max] = size(X);
r = size(U_tr,2);
T = Indicator_full;
PER.EV    = zeros(1,t_max);
PER.SEP   = zeros(1,t_max);
PER.SE    = zeros(1,t_max);
PER.Angle = zeros(1,t_max);

r_0 = r;
ev_thresh = 7.5961e-04;
K     = 32;
alpha = 60;
R     = 0;
overlap_step = alpha;
hop          = overlap_step; 

P_init = orth(randn(n, r));
P_hat  = P_init;

S_hat        = zeros(n, t_max);
S_hat_buffer = zeros(n, t_max);
L_hat        = zeros(n, t_max);
L_hat_buffer = zeros(n, t_max);
t_hat        = [];

k   = 0;
cnt = 1;
P_track_full{cnt} = P_hat;
T_calc(cnt)       = 1;
ph = 1;     %ph - 0 => detect, 1 => pca
nextbatch_start = 0;
alpha0 = alpha;
permit = 0;

for ii = 1: t_max     
    
    %% Estimate support
    Atf.times = @(x) x - (P_hat * (P_hat' * x));
    Atf.trans = @(y) y - (P_hat * (P_hat' * y));
    phi.times = @(x) x - (P_hat * (P_hat' * x));
    y_t = Atf.times(X(:, ii));
    if (ii == 1)
        x0 = zeros(n,1);
    else
%         x0 = S_hat(:,ii-1);
    end
    
    T_union = find(T(:, ii) == 0);
    tol = 1e-16;
    S_hat(T_union, ii) = ccgls(@Phifun, y_t, T_union, P_hat, ...
        0, tol, 20);
    L_hat(:, ii) = X(:, ii) - S_hat(:, ii);    
    
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

                for reuse = 1:R
                    for kk = idx_buffer
                        Atf.times = @(x) x - (P_hat * (P_hat' * x));
                        Atf.trans = @(y) y - (P_hat * (P_hat' * y));
                        y_t_buffer = Atf.times(X(:, kk));
                        

                        T_union_buffer = find(T(:,kk) == 0);
                        S_hat_buffer(T_union_buffer, kk) = ccgls(@Phifun, y_t_buffer, ...
                            T_union_buffer, P_hat, ...
                            0, tol, 20);
                        L_hat_buffer(:, kk) = X(:, kk) - S_hat_buffer(:, kk);                        
                    end
                    P_hat = simpleEVD(L_hat_buffer(:, idx_buffer), r);
                end
            end
                        
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
    PER.SEP(1,ii)   = sub_est_per(U_tr,P_hat,'SEP');
    PER.Angle(1,ii) = sub_est_per(U_tr,P_hat,'Angle');
    PER.SE(1,ii)    = sub_est_per(U_tr,P_hat,'SE');
    PER.EV(1,ii)    = sub_est_per(U_tr,P_hat,'EV'); 

end
U_est = P_hat; 
end