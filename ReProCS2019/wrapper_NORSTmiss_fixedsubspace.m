clear;
clc;

addpath('PROPACK')

%% Algorithms to run
NORST = 1;
NORST_OFFLINE = 0;

%% Parameter Initialization
n = 1000;
t_max = 3000;
alpha = 60;
f = 100;
MC = 1;
t_calc_pca = [1:alpha: t_max];
t_calc = t_calc_pca;

%NORST
temp_SE_NORST = zeros(length(t_calc_pca), MC);
temp_err_L_NORST = zeros(t_max, MC);
t_NORST = 0;
err_L_fro_NORST = zeros(MC,1);

temp_SE_NORST_off = zeros(length(t_calc_pca), MC);
temp_err_L_NORST_off = zeros(t_max, MC);
t_NORST_off = 0;
err_L_fro_NORST_offline = zeros(MC,1);

for mc = 1 : MC
    fprintf('Monte-Carlo iteration %d in progress \n', mc);
    
    % b. Bernoulli Model
    rho = 0.1; % fraction of missing entries
    BernMat = rand(n, t_max);
    T = 1 .* (BernMat <= 1 - rho); % observed entries' support
    
    %% Generating low-rank matrix
    r_0 = 30;
    L = zeros(n, t_max);
    
    lambda_min = sqrt(f)/2;
    lambda_max = sqrt(f);
    
    offset = 0; %if offset is not zero, eigenvalues (Lambda) are varying in time;
    diag_entries1 = offset +  [linspace(lambda_max, lambda_min, r_0)];
    diag_entries2 = -offset + [linspace(lambda_max, lambda_min, r_0)];
    
    coeff_train = zeros(r_0, t_max);
    for cc = 1 : r_0
        coeff_train(cc, 1:2:end-1) = -diag_entries1(cc) + ...
            2 * diag_entries1(cc) * rand(1, t_max/2);
        
        coeff_train(cc, 2:2:end) = -diag_entries2(cc) + ...
            2 * diag_entries2(cc) * rand(1, t_max/2);
    end
    
    P = orth(randn(n, r_0));
    L = P * coeff_train;
    
    eps_noise = 0; % noise
    L = L + eps_noise * (rand(n,t_max) - 0.5);
    M = L .* T ;
    
    %% Algorithm parameters for NORST
    if(NORST == 1)
        fprintf('\tNORST\t');
        K = 32;
        ev_thresh = 7.5961e-04;
        tol = 1e-16;
        overlap_step = alpha; % if it is set to alpha then windows don't overlap
        R = 0; % number of reuse
        
        %     P_init = orth(randn(n,r_0));
        P_init = zeros(n,r_0);
        
        t_norst = tic;
        [L_hat, P_hat, S_hat, t_hat, P_track_full, t_calc] =  ...
            NORST_random(M, T, r_0, ev_thresh, alpha,K,R,overlap_step);
        t_NORST = toc(t_norst)
        err_L_fro_NORST(mc) = norm(L-L_hat,'fro')/norm(L,'fro');
        
    end
    if(NORST_OFFLINE == 1)
        fprintf('\tNORST-offline\t');
        ev_thresh = 7.5961e-04;
        tol = 1e-16;
        K_off = 32;
        [BG,L_hat_off, P_hat_off, S_hat_off, t_hat_off, P_track_full_off, t_calc_off] =  ...
            NORST_offline(M, T, r_0, tol, ev_thresh, alpha, K_off);
        %     [BG_off, L_hat_off, P_hat_off, S_hat_off, t_hat_off, P_track_full_off, t_calc_off] =  ...
        %         NORST_random_offline(M, T, r_0, tol, ev_thresh, alpha, K_off, buffer);
        err_L_fro_NORST_offline(mc) = norm(L-L_hat_off,'fro')/norm(L,'fro');
    end
    
    %% Compute Performance Metrics
    %compute the "frobenius norm errors"
    if(NORST == 1)
        temp_err_L_NORST(:, mc) = sqrt(mean((L - L_hat).^2, 1)) ...
            ./ sqrt(mean(L.^2, 1));
    end
    if(NORST_OFFLINE == 1)
        temp_err_L_NORST_off(:, mc) = sqrt(mean((L - L_hat_off).^2, 1)) ...
            ./ sqrt(mean(L.^2, 1));
    end
    
    %computing subspace errors
    for jj = 1 : length(t_calc_pca)
        if(NORST == 1)
            temp_SE_NORST(jj, mc) = ...
                Calc_SubspaceError(P_track_full{jj}, P);
        end
        
        if(NORST_OFFLINE == 1)
            temp_SE_NORST_off(jj, mc) = ...
                Calc_SubspaceError(P_track_full_off{1}, P);
        end
    end
    fprintf('\n')
end

err_SE_NORST = mean(temp_SE_NORST, 2);
err_L_NORST = mean(temp_err_L_NORST, 2);

err_SE_NORST_off = mean(temp_SE_NORST_off, 2);
err_L_NORST_off = mean(temp_err_L_NORST_off, 2);

figure
strx = 't';
stry = '$$\log_{10} dist(\hat{P}_{(t)}, P_{(t)})$$';
%
semilogy(t_calc_pca(1:2:end),err_SE_NORST(1:2:end),'-*r','LineWidth',2,'MarkerSize',10);
grid on

xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 20);