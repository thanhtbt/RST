%%%Wrapper to obtain a phase transition matrix to show the "constant"
%%%fraction of outlier bounds for NORST while showing standard RPCA
%%%techniues can only tolerate 1/r.

%clear;
%clc;
close all

addpath('YALL1_v1.4')

%% Parameter Initialization
n = 1000;
t_max = 11000;
t_train = 100;
alpha = 300;
f = 100;
MC = 100;

t_calc_pca = [alpha - 1 : alpha : t_max - t_train];

outfracrowrange = linspace(0.01, 0.7, 10);
rrange = ceil(linspace(3, 30, 10));

PhaseTransNORST = zeros(length(outfracrowrange), length(rrange), MC);
PhaseTransAltProj = zeros(length(outfracrowrange), length(rrange), MC);

ttall = tic;

for bb  = 1 : length(outfracrowrange)
    b0 = outfracrowrange(bb);
    for rr = 1 : length(rrange)
        r_0 = rrange(rr);
        parfor mc = 1 : MC
            fprintf('Monte-Carlo iteration %d\t b0 %.2f\t r %d', mc, b0, r_0);
            
            %%%Generating support set and sparse vectors
            %Using bernoulli model since AltProj fails on Moving Object
            T = zeros(n, t_max);
            rho_train = 0.02;
            rho = b0;
            x_min = 10;
            x_max = 20;
            
            BernMat = rand(n, t_max);
            T(:, 1 : t_train) = 1.* (BernMat(:, 1 : t_train) <= rho_train);
            T(:, t_train + 1 : end) = 1 .* (BernMat(:, t_train + 1 : t_max) <= rho);
            S = (x_min + (x_max - x_min) * rand(n, t_max)) .* T;
            
            %%%Generate low-rank matrix
            %r_0 = 30;
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
            
            %%Generate Subspaces
            Btemp1 = randn(n);
            B1 = (Btemp1 - Btemp1')/2;
            Btemp2 = randn(n);
            B2 = (Btemp2 - Btemp2')/2;
            
            delta1 = .5e-2;
            delta2 = 0.8 * delta1;
            P = orth(randn(n, r_0));
            PP1 = expm(delta1 * B1)  * P;
            PP2 = expm(delta2 * B2) * PP1;
            
            L(:, 1:t_1) = P(:, 1:r_0) * coeff_train(:, 1:t_1);
            L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
            L(:, t_2 + 1 : end) = PP2 * coeff_train(:, t_2+1:end);
            M = L + S;
            
            %% Calls to RPCA algorithms
            
            %%%Algorithm parameters for NORST
            K = 8;
            omega = x_min / 2;
            ev_thresh = 7.5961e-04;
            
            %%%Call to offline NORST
            fprintf('\tOffline NORST\t');
            P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 15));
            [L_hat_off, P_hat_off, S_hat_off, T_hat_off, t_hat_off, ...
                P_track_full_off, P_track_new_off] = Offline_NORST(...
                M(:, t_train + 1 :end), P_init, ev_thresh, alpha, K, omega);
            
            %%Call to AltProj
            fprintf('AltProj\n');
            [L_hat_ncrpca, S_hat_ncrpca] = ...
                ncrpca(M(:, t_train + 1 : end), 3 * r_0, 1e-6, 100);
            
            %% Compute performance metrics
%             temp_err_L_ncrpca = sqrt(mean((L(:, t_train + 1 :end) - ...
%                 L_hat_ncrpca).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
%             temp_err_L_norst = sqrt(mean((L(:, t_train + 1 :end) - ...
%                 L_hat_off).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
            
            PhaseTransAltProj(rr, bb, mc) = ...
                norm(L(:, t_train + 1 : end) - L_hat_ncrpca, 'fro') / ...
                norm(L(:, t_train + 1 : end), 'fro');
            PhaseTransNORST(rr, bb, mc) = ...
                norm(L(:, t_train + 1 : end) - L_hat_off, 'fro') / ...
                norm(L(:, t_train + 1 : end), 'fro');
            
%             PhaseTransNORST(rr, bb, mc) = mean(temp_err_L_norst);
%             PhaseTransAltProj(rr, bb, mc) = mean(temp_err_L_ncrpca);
        end
    end
end
fprintf('\n')
toc(ttall)

save('phase_trans_b0_r_LF.mat');

% save('data_TIT/phase_trans_b0_r.mat');