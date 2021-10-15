clear;
clc;

addpath('PROPACK')

%% Algorithms to run
NORST = 1;

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

for mc = 1 : MC
     fprintf('Monte-Carlo iteration %d in progress \n', mc);

    %% Generating support set using: a. Moving Object Model    b. Bernoulli Model
    % a. Moving Object Model
%     s = 100;            
%     S = zeros(n, t_max);
%     rho_s = 1;
%     b0 = 0.1;
%     beta = ceil(b0 * alpha1);
%     x_max = 25;
%     x_min = 15;
%     alpha1 = 100;
%     num_changes = floor((t_max -t_train)/beta);
% 
%     num_changes1 = min(floor(alpha1 / beta), ceil(n/s));
% 
%     flag = 0;
%     ii1 = 1;
%     fval1 = 0;
%     for ii = 1 : num_changes
%         if(~flag)   %%downward motion
%             if(ii1 <= num_changes1)
%                 bind = fval1 + (ii1 - 1) * s/rho_s + 1;
%                 sind = min(bind - 1 + s, n);
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 1;
%                     ii1 = 1;
%                     fval2 = bind;
%                 end
%             end
%         else
%             if(ii1 <= num_changes1)
%                 bind = max(fval2 - (ii1 - 1) * s/rho_s , 1);
%                 sind = bind - 1 + s;
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 0;
%                     ii1 = 1;
%                 end
%             end
%         end
%         idx = bind : sind;
%         jdx = (ii-1) * beta + 1 : ii * beta;
%         S(idx, jdx) = x_min + ...
%             (x_max - x_min) * rand(length(idx), beta);
%         T(idx, jdx) = 1;
%     end
%     
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
    delta_t = 100;
    U0 = P;
    subspace_size = 1500;
    U_track = cell(ceil(t_max/subspace_size),1);
    for i=1:length(U_track)
        Btemp1 = randn(n);
        B1 = (Btemp1 - Btemp1')/2;
        
        t_1 = (i-1)*subspace_size + 1;
        t_2 = min(i*subspace_size,t_max);
        
        U_track{i} = U0;
        L(:, t_1:t_2) = U0 * coeff_train(:,t_1:t_2);
        
        % subspace change at each time t
%         if t == 1499
%             U = expm(100*B1)*U0;        
%         else
%             delta_t = 1e-7;
%             U = expm(delta_t*B1)*U0;
%         end 
        
        U = expm(delta_t*B1)*U0;        
        U0 = U;
    end 
    eps_noise = 1e-3; % noise
    L = L + eps_noise * (rand(n,t_max) - 0.5);
    M = L .* T ;
    
    %% Algorithm parameters for NORST
    if(NORST == 1)
    fprintf('\tNORST\t');
    K = 25;
    ev_thresh = 7.5961e-04;
    tol = 1e-16;
    overlap_step = alpha; % if it is set to alpha then windows don't overlap
    R = 0; % number of reuse 
    
%     P_init = orth(randn(n,r_0));
    P_init = zeros(n,r_0);
    
    t_norst = tic;
    [L_hat, P_hat, S_hat, t_hat, P_track_full, t_calc] =  ...
        NORST_random(M, T, r_0, ev_thresh, alpha, K,R,overlap_step);
    t_NORST = toc(t_norst)
    err_L_fro_NORST(mc) = norm(L-L_hat,'fro')/norm(L,'fro');
    
    end

    %% Compute Performance Metrics
    %compute the "frobenius norm errors"
if(NORST == 1)
    temp_err_L_NORST(:, mc) = sqrt(mean((L - L_hat).^2, 1)) ...
        ./ sqrt(mean(L.^2, 1));
end

    %computing subspace errors
   for jj = 1 : length(t_calc_pca)
            tt = ceil(t_calc_pca(jj)/subspace_size);
            if(NORST == 1)           
            temp_SE_NORST(jj, mc) = ...
                Calc_SubspaceError(P_track_full{jj}, U_track{tt});
            end
            
            if(NORST_OFFLINE == 1)           
            temp_SE_NORST_off(jj, mc) = ...
                Calc_SubspaceError(P_track_full{jj}, U_track{tt});
            end
    end
fprintf('\n')
end

err_SE_NORST = mean(temp_SE_NORST, 2);
err_L_NORST = mean(temp_err_L_NORST, 2);

figure
strx = 't';
stry = '$$\log_{10} dist(\hat{P}_{(t)}, P_{(t)})$$';
% 
semilogy(t_calc_pca(1:2:end),err_SE_NORST(1:2:end),'-*r','LineWidth',2,'MarkerSize',10);
grid on

xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 20);