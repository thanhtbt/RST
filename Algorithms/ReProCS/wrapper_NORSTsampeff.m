clear;
clc;

addpath('PROPACK')

%% Algorithms to run
NORST = 1;
NORSTeff = 1;

%% Parameter Initialization
n = 1000;
t_max = 3000;
alpha = 60;
f = 100;
MC = 1;
t_calc_pca = [1:alpha: t_max];

%NORST
temp_SE_NORST = zeros(length(t_calc_pca), MC);
temp_err_L_NORST = zeros(t_max, MC);
t_NORST = 0;
err_L_fro_NORST = zeros(MC,1);

%NORSTeff
temp_SE_NORSTeff = zeros(length(t_calc_pca), MC);
temp_err_L_NORSTeff = zeros(t_max, MC);
t_NORSTeff = 0;
err_L_fro_NORSTeff = zeros(MC,1);

for mc = 1 : MC
%     fprintf('(rho2 = %0.1f) Monte-Carlo iteration %d\t',rho2,mc);  
    %%% bernoulli model for observed entries
            rho1 = 0.1; %denotes fraction of number of missing entries
            rho2 = 0.85;            
            Ka = 900;
            
            BernMat1 = rand(n, Ka);
            BernMat2 = rand(n, t_max - Ka);
            T1 = 1 .* (BernMat1 <= 1 - rho1);
            T2 = 1 .* (BernMat2 <= 1 - rho2);
            T = [T1,T2];
    
    %%%Generate low-rank matrix
    r_0 = 30;
    
    lambda_min = sqrt(f)/2;
    lambda_max = sqrt(f);
    
    offset = 0;
    
    diag_entries1 = offset +  [linspace(lambda_max, lambda_min, r_0)];
    diag_entries2 = -offset + [linspace(lambda_max, lambda_min, r_0)];
    t_1 = t_max;
    
    coeff_train = zeros(r_0, t_max);
    for cc = 1 : r_0
        coeff_train(cc, 1:2:end-1) = -diag_entries1(cc) + ...
            2 * diag_entries1(cc) * rand(1, t_max/2);
        
        coeff_train(cc, 2:2:end) = -diag_entries2(cc) + ...
            2 * diag_entries2(cc) * rand(1, t_max/2);
    end
    
    %%Generate the Subspace

    P = orth(randn(n, r_0));
    L = P * coeff_train(:, 1:t_1);

    M = L .* T ;
    
    %% Calling MC/ST algorithms
    %%Algorithm parameters for NORST
    if(NORST == 1)
    fprintf('\tNORST\t');
    K = 33;
    ev_thresh = 7.5961e-04;
    tol = 1e-16;
    overlap_step = alpha;
    R = 0;
%     P_init = orth(randn(n,r_0));
    P_init = zeros(n,r_0);
    t_norst = tic;
    [L_hat, P_hat, S_hat, t_hat, P_track_full, t_calc] =  ...
        NORST_random(M, T, r_0, ev_thresh, alpha, K,R,overlap_step);
        
    L_hat_norst = L_hat;

    t_NORST = toc(t_norst);
  
    err_L_fro_NORST(mc) = norm(L-L_hat_norst,'fro')/norm(L,'fro');
    end
    
    if(NORSTeff == 1)
    fprintf('\tNORST-samp-eff\t');
    K = 15;
    ev_thresh = 7.5961e-04;
    tol = 1e-16;
    overlap_step = alpha;
    R = 0;
%     P_init = orth(randn(n,r_0));
    P_init = zeros(n,r_0);
    t_norsteff = tic;
    [L_hat, P_hat, S_hat, t_hat, P_track_full, t_calc] =  ...
        NORST_random(M(:,1:Ka), T1, r_0, ev_thresh, alpha, K,R,overlap_step);
    
    L_hat_eff = NORSTefficient(P_hat,M(:,Ka+1:end),T2);
    L_hat_fin = [L_hat,L_hat_eff];

    t_NORSTeff = toc(t_norsteff);
  
    err_L_fro_NORSTeff(mc) = norm(L-L_hat_fin,'fro')/norm(L,'fro');
    end


    %% Compute Performance Metrics
    %compute the "frobenius norm errors"
if(NORST == 1)
    temp_err_L_NORST(:, mc) = sqrt(mean((L - L_hat_norst).^2, 1)) ...
        ./ sqrt(mean(L.^2, 1));
end

if(NORSTeff == 1)
    temp_err_L_NORSTeff(:, mc) = sqrt(mean((L - L_hat_fin).^2, 1)) ...
        ./ sqrt(mean(L.^2, 1));
end

fprintf('\n')
end

err_L_NORST = mean(temp_err_L_NORST, 2);
err_L_NORSTeff = mean(temp_err_L_NORSTeff, 2);

figure
strx = 't';
stry = '$$\log_{10} dist(\hat{P}_{(t)}, P_{(t)})$$';
 
p1 = semilogy(1:alpha:t_max,err_L_NORST(1:alpha:t_max),'-*r','LineWidth',2,'MarkerSize',10);
hold on
grid on
p2 = semilogy(1:alpha:t_max,err_L_NORSTeff(1:alpha:t_max),'-sg','LineWidth',2,'MarkerSize',10);
xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 20);
legend([p1,p2],{'NORST-miss','NORST-samp-eff'},'FontSize',20)