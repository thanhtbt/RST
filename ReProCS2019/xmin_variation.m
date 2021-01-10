%%%Wrapper to verify the variation of NORST and ALtProj with x_min, the
%%%smallest outlier-magnitude.

clear;
clc;
close all

addpath('YALL1_v1.4')

%% Parameter Initialization
n = 1000;
t_max = 12000;
t_train = 100;
miss_s = 0;
miss_s_pca = 0;
alpha = 300;
f = 100;
MC = 3;

xminBnd = [.5, 5, 10];
temp_SE_reprocs_pca = zeros(length(xminBnd), ceil((t_max- t_train)/alpha) - 1, MC);
temp_SE_reprocs_off = zeros(length(xminBnd), ceil((t_max- t_train)/alpha) - 1, MC);
temp_SE_ncrpca = zeros(length(xminBnd), ceil((t_max- t_train)/alpha) - 1, MC);


ttall = tic;
for mc = 1 : MC
    temp_err_L_pca = zeros(length(xminBnd), t_max - t_train);
    temp_err_L_off = zeros(length(xminBnd), t_max - t_train);
    temp_err_L_ncrpca = zeros(length(xminBnd), t_max- t_train);
    for xx = 1 : length(xminBnd)
        x_min = xminBnd(xx);
        
        fprintf('Monte-Carlo iteration %d in progress for xmin %.2f\n', mc, x_min);
        
        
        %%%Generating support set and sparse vectors
        S = zeros(n, t_max);
        T = zeros(n, t_max);
        
        b0 = 0.35;
        rho = 1;
        alpha1 = 100;
        s = 50;
        s_train = s/5;
        alpha1_train = alpha1/2;
        beta = ceil(b0 * alpha1);
        x_max = 20;
        %x_min = 10;
        alpha_train = alpha1;
        num_changes = floor((t_max - t_train)/beta);
        
        num_changes1 = min(floor(alpha1 / beta), ceil(n/s));
        
        %training small outlier fraction
        beta_train = 1;
        num_train = floor(t_train/ beta_train);
        num_train_1 = min(floor(alpha1_train / beta_train), ceil(n/s_train));
        fval1 = 0;
        flag = 0;
        ii1 = 1;
        for ii = 1 : num_train
            if(~flag) %%downward motion
                if(ii1 <= num_train_1)
                    bind = fval1 + (ii1 - 1) * s_train/rho + 1;
                    sind = min(bind - 1 + s_train, n);
                    ii1 = ii1 + 1;
                    if(ii1 == num_train_1 + 1)
                        flag = 1;
                        ii1 = 1;
                        fval2 = bind;
                    end
                end
            else
                if(ii1 <= num_train_1)
                    bind = max(fval2 - (ii1 - 1) * s_train/rho, 1);
                    sind = bind - 1 + s_train;
                    ii1 = ii1 + 1;
                    if(ii1 == num_train_1 + 1)
                        flag = 0;
                        ii1 = 1;
                    end
                end
            end
            idx = bind : sind;
            jdx = (ii-1) * beta_train + 1 : ii * beta_train;
            S(idx, jdx) = x_min + ...
                (x_max - x_min) * rand(length(idx), beta_train);
            T(idx, jdx) = 1;
        end
        
        flag = 0;
        ii1 = 1;
        fval1 = 0;
        
        for ii = 1 : num_changes
            if(~flag)   %%downward motion
                if(ii1 <= num_changes1)
                    bind = fval1 + (ii1 - 1) * s/rho + 1;
                    sind = min(bind - 1 + s, n);
                    ii1 = ii1 + 1;
                    if(ii1 == num_changes1 + 1)
                        flag = 1;
                        ii1 = 1;
                        fval2 = bind;
                    end
                end
            else
                if(ii1 <= num_changes1)
                    bind = max(fval2 - (ii1 - 1) * s/rho , 1);
                    sind = bind - 1 + s;
                    ii1 = ii1 + 1;
                    if(ii1 == num_changes1 + 1)
                        flag = 0;
                        ii1 = 1;
                    end
                end
            end
            idx = bind : sind;
            jdx = t_train + (ii-1) * beta + 1 : t_train + ii * beta;
            %setting all entries to x_min to make things harder!
            S(idx, jdx) = x_min * ones(length(idx), beta);
            T(idx, jdx) = 1;
        end
        
        %%%Generate low-rank matrix
        r_0 = 30;
        r_1 = 0;
        r_2 = 0;
        r = r_0;
        L = zeros(n, t_max);
        
        diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0)];
        t_1 = 3000;
        t_2 = 8000;
        coeff_train = zeros(r_0, t_max);
        
        for cc = 1 : r_0
            coeff_train(cc, :) = -diag_entries(cc) + ...
                2 * diag_entries(cc) * rand(1, t_max);
        end
        
        Btemp1 = randn(n);
        B1 = (Btemp1 - Btemp1')/2;
        Btemp2 = randn(n);
        B2 = (Btemp2 - Btemp2')/2;
        
        delta1 = .5e-3;
        delta2 = 0.8 * delta1;
        
        P = orth(randn(n, r_0));
        PP1 = expm(delta1 * B1)  * P;
        PP2 = expm(delta2 * B2) * PP1;
        
        L(:, 1:t_1) = P(:, 1:r_0) * coeff_train(:, 1:t_1);
        L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
        L(:, t_2 + 1 : end) = PP2 * coeff_train(:, t_2+1:end);
        M = L + S;
        
        %% Calls to NORST, Offline NORST and AltProj
        
        %%%Algorithm parameters
        K = 8;
        omega = x_min / 2;
        %     gamma = sqrt(4 * log(n)/n);
        %     s = ceil((gamma + rho) * n);
        
        %%%Call to NORST
        fprintf('NORST\t');
        ev_thresh = 7.5961e-04;
        P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 15));
        [L_hat_pca, P_hat_pca, S_hat_pca, T_hat_pca, t_hat_pca, ...
            P_track_full_pca, t_calc_pca] = ...
            NORST(M(:, t_train + 1 :end),...
            P_init, ev_thresh, alpha, K, omega);
        
        %%%Call to offline NORST
        fprintf('Offline NORST\n');
        P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 15));
        [L_hat_off, P_hat_off, S_hat_off, T_hat_off, t_hat_off, ...
            P_track_full_off, P_track_new_off] = ...
            Offline_NORST(M(:, t_train + 1 :end),...
            P_init, ev_thresh, alpha, K, omega);
        
        %%Call to AltProj
        fprintf('AltProj\t');
        [L_hat_ncrpca, S_hat_ncrpca] = ...
            ncrpca(M(:, t_train + 1 : end), r_0 + r_1 + r_2, 1e-6, 100);
        
        %%Compute performance metrics
        temp_err_L_pca(xx, :) = ...
            sqrt(mean((L(:, t_train+1:end) - L_hat_pca).^2, 1)) ...
            ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
        temp_err_L_off(xx, :) = ...
            sqrt(mean((L(:, t_train+1:end) - L_hat_off).^2, 1)) ...
            ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
        temp_err_L_ncrpca(xx, :) = sqrt(mean((L(:, t_train + 1 :end) - ...
            L_hat_ncrpca).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
        
        %     miss_s = ...
        %         miss_s + (length(find(S_hat_off))- length(find(S)))/numel(S);
        %     miss_s_pca = ...
        %         miss_s_pca + (length(find(S_hat_pca))- length(find(S)))/numel(S);
        
        
        
        for jj = 1 : length(t_calc_pca)
            P_hat_ncrpca = orth(L_hat_ncrpca(:, 1 : t_calc_pca(jj)));
            if (t_calc_pca(jj) +t_train < t_1)
                temp_SE_reprocs_pca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_pca{jj}, ...
                    P(:, 1:r_0));
                temp_SE_reprocs_off(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_off{1}, ...
                    P(:, 1:r_0));
                temp_SE_ncrpca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_hat_ncrpca, ...
                    P(:, 1 :r_0));
            elseif((t_calc_pca(jj) +t_train >= t_1) && (t_calc_pca(jj) + t_train < t_2))
                temp_SE_reprocs_pca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_pca{jj}, PP1);
                temp_SE_reprocs_off(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_off{2}, PP1);
                temp_SE_ncrpca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_hat_ncrpca, PP1);
            else
                temp_SE_reprocs_pca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_pca{jj}, PP2);
                temp_SE_reprocs_off(xx, jj, mc) = ...
                    Calc_SubspaceError(P_track_full_off{3}, PP2);
                temp_SE_ncrpca(xx, jj, mc) = ...
                    Calc_SubspaceError(P_hat_ncrpca, PP2);
            end
        end
        
    end
end
fprintf('\n')
toc(ttall)

save('data_TIT/xmin_var_mc3_xmaxhigh.mat')