clear; clc; 

addpath('YALL1_v1.4')

nb_exp = 5;
outlier_fac = 1;
outlier_spr =  0.05;
SAMPLING    =  0.9;
std_s     = 1;
fac_noise = 1;
snr       = 20;
std_brt   = fac_noise * 10^(-snr/20);

numr = 50; numc = 5000;
r = 2;

S       = std_s * randn(r,numc);
A       = randn(numr,r);
X0      = A*S;
V       = orth(X0);
V       = V(:,1:r);

normX = norm(X0,'fro');


% Missing value
ind = randsrc(numr,numc,[0 1; (1-SAMPLING) SAMPLING]); ind_vec = ind(:)';
I = [];
J = [];
for ii = 1 : numc
    locs = find(ind(:,ii)==1);
    I = [I locs'];
    jj = ii*ones(1,length(locs));
    J = [J jj];
end
eta_1  = zeros(1,numc);      rho_1 = zeros(1,numc);
eta_2  = zeros(1,numc);      rho_2 = zeros(1,numc);
eta_3  = zeros(1,numc);      rho_3 = zeros(1,numc);
eta_4  = zeros(1,numc);      rho_4 = zeros(1,numc);
eta_5  = zeros(1,numc);      rho_5 = zeros(1,numc);
eta_6  = zeros(1,numc);      rho_6 = zeros(1,numc);
eta_7  = zeros(1,numc);      rho_7 = zeros(1,numc);
eta_8  = zeros(1,numc);      rho_8 = zeros(1,numc);
eta_9  = zeros(1,numc);      rho_9 = zeros(1,numc);
eta_10  = zeros(1,numc);     rho_10 = zeros(1,numc);

% ============================== PROGRAMME  ==============================

for ii_exp = 1: nb_exp
    fprintf('\nExperience %d of %d ===========================================\n', ii_exp,nb_exp)
    
    
    N   =  randn(numr, numc);
    normN = norm(N,'fro');
    X   = X0 + std_brt * normX/normN * N;
    
    X_free = X; %free outliers
    
    C   = zeros(numr,numc);
    p   = randperm(numr*numc);
    L   = round(outlier_spr*numr*numc);
    mgB = max(abs(X(:)));
    C(p(1:L))  = outlier_fac*mgB * randn(L,1);
    
    X = X + C;
    X_free_missing      = X;
    XzeroOutliers       = X;
    XzeroOutliers(C~=0) = 0;
    

    %% Algorithms
    
    fprintf('+ NORST: ')
    t_start = tic;
    OPTS  = [];
    [rho, eta]   = AutoReProCS_ST(X,V,OPTS);
    toc(t_start)
    rho_10 = rho_10 + rho;
    eta_10 = eta_10 + eta;
    % figure ;  semilogy(rho)
    
end
k = 1;
figure; semilogy(rho_10(1:k:numc)/nb_exp)





%%




