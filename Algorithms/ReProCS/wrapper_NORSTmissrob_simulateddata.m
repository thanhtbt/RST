clear;
clc;

addpath('YALL1_v1.4')
addpath('PROPACK')

n = 1000; % dimensionality
t_max = 10000; % number of time instances
alpha = 60; % frame length
f = 100; % condition number

t_calc_pca = alpha-1:alpha:t_max;
temp_SE_NORST = zeros(length(t_calc_pca), MC);
temp_err_L_NORST = zeros(t_max, MC);
err_L_fro = zeros(MC,1);

cnt = 1;
MC = 100;
for mc = 1 : MC            
    fprintf('Monte-Carlo iteration %d in progress \n', mc);

    %% Generating support set and sparse vectors
    s = 100;            
    S = zeros(n, t_max);
    rho_s = 1;
    b0 = 0.1;
    beta = ceil(b0 * alpha1);
    x_max = 25;
    x_min = 15;
    alpha1 = 100;
    num_changes = floor((t_max -t_train)/beta);

    num_changes1 = min(floor(alpha1 / beta), ceil(n/s));

    flag = 0;
    ii1 = 1;
    fval1 = 0;
    for ii = 1 : num_changes
        if(~flag)   %%downward motion
            if(ii1 <= num_changes1)
                bind = fval1 + (ii1 - 1) * s/rho_s + 1;
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
                bind = max(fval2 - (ii1 - 1) * s/rho_s , 1);
                sind = bind - 1 + s;
                ii1 = ii1 + 1;
                if(ii1 == num_changes1 + 1)
                    flag = 0;
                    ii1 = 1;
                end
            end
        end
        idx = bind : sind;
        jdx = (ii-1) * beta + 1 : ii * beta;
        S(idx, jdx) = x_min + ...
            (x_max - x_min) * rand(length(idx), beta);
        T(idx, jdx) = 1;
    end
    
    fprintf('fraction of sparse entries: %d \n',length(find(T(:) == 1)) / numel(T));
    
    t_train = 400;
    
    %% Missing entries' support
    % Bernoulli Model
    rho = 0.1; % fraction of missing entries
    BernMat = rand(n, t_max);
    T_obs = 1 .* (BernMat <= 1 - rho); % observed entries' support 
    
    %% Generating low-rank matrix
    L = zeros(n, t_max);
    r = 30;
    diag_entries = linspace(sqrt(f), sqrt(f)/2, r);
    
    % direction change times
    t_1 = 4000;
    t_2 = 8000;
    
    coeff_train = zeros(r, t_max);
    for cc = 1 : r
        coeff_train(cc, :) = -diag_entries(cc) + ...
            2 * diag_entries(cc) * rand(1, t_max);
    end

    Btemp1 = randn(n);
    B1 = (Btemp1 - Btemp1')/2;
    Btemp2 = randn(n);
    B2 = (Btemp2 - Btemp2')/2;

    delta1 = 0.5e-3;
    delta2 = 0.8 * delta1;

    P = orth(randn(n, r));
    PP1 = expm(delta1 * B1)  * P;
    PP2 = expm(delta2 * B2) * PP1;

    L(:, 1:t_1) = P(:, 1:r) * coeff_train(:, 1:t_1);
    L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
    L(:, t_2 + 1 : end) = PP2 * coeff_train(:, t_2+1:end);
    
    M = L + S;
    M = M .* T_obs;    
    M_norst = M(:,t_train+1:end);
    
    %% Calls to NORST
    fprintf('NORST-miss-robust\n');
    % Algorithm parameters
    K = 36;
    omega = x_min / 2;   
    ev_thresh = 7.5961e-04;
    tol = 1e-16;
    
    M_train = M(:,1:t_train);
    M_train(M_train == 0)= x_max;
    
    t_norst = tic;
    P_init = orth(ncrpca(M_train, r, 1e-3, 100));
    
    [L_hat, P_hat, S_hat, T_hat,...
     t_hat,P_track_full, t_calc_norst] = NORST(M_norst,...
               T_obs(:, 1+t_train:end),P_init,ev_thresh,alpha,K,omega,tol);
    
    t_NORST = toc(t_norst);
    
    
    
    %% Compute performance metrics
    % mean least square error in the low-rank matrix
    temp_err_L(mc, :) = ...
        sqrt(mean((L(:, t_train + 1 : end) - L_hat).^2, 1)) ./ ...
        sqrt(mean(L(:, t_train + 1 : end).^2, 1));
    % frobenious norm
    err_L_fro(mc) = norm(L(:,t_train+1:end)-L_hat,'fro')/norm(L(:,t_train+1:end),'fro');        
        
    % subspace error
    for jj = 1 : length(t_calc_norst)
        if(t_calc_norst(jj) < t_1)            
            temp_SE_Phat_P(mc, jj) = ...
                Calc_SubspaceError(P_track_full{jj}, P);
        
        elseif((t_calc_norst(jj) >= t_1) && (t_calc_norst(jj) < t_2))            
            temp_SE_Phat_P(mc, jj) = ...
                Calc_SubspaceError(P_track_full{jj}, PP1);            
        else            
            temp_SE_Phat_P(mc, jj) = ...
                Calc_SubspaceError(P_track_full{jj}, PP2);
        end
    end
end

err_L = mean(temp_err_L, 1);
SE_Phat_P = mean(temp_SE_Phat_P, 1);
 
str1 = 't';
str2 = '$$\log SE(\hat{P}, P)$$';
str3 = ['\rho_obs = ',num2str(1-rho)];
str4 = '$$\log \frac{||\hat{l}-l||^2}{||l||^2}$$';

figure
plot(t_calc_norst + t_train,log10(SE_Phat_P_rpca),'r*--','LineWidth',1,'MarkerSize',6)
grid on
xlabel(str1,'interpreter', 'latex','FontSize',20)
ylabel(str2,'interpreter', 'latex','FontSize',20)
title(str3);

figure
step_size = 2*alpha;
plot(t_train + 1 :step_size: t_max, log10(err_L(1:step_size:t_max-t_train)),'b*--','LineWidth',1,'MarkerSize',6)
xlabel(str1,'interpreter', 'latex','FontSize',20)
ylabel(str4,'interpreter', 'latex','FontSize',20)
