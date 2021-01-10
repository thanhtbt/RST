clear; clc;

close all; 



nb_exp = 1;
outlier_fac = 1;
outlier_den =  0.002;
SAMPLING    =  0.2;
std_s     = 1;
fac_noise = 1;
snr     = 20;
std_brt = fac_noise * 10^(-snr/20);

numr = 50; numc = 5000;
r = 5;

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

for ii_exp = 1: nb_exp
    fprintf('\nExperience %d of %d ===========================================\n', ii_exp,nb_exp)

    S       = std_s * randn(r,numc);
    A       = randn(numr,r);
    X0      = A*S; % X0 = X0 / norm(X0);
    V       = orth(X0);
    normX   = norm(X0,'fro'); 
    N       = randn(numr, numc); 
    normN = norm(N,'fro');
    X       = X0 + std_brt * (normX/normN) * N;
    X_free  = X;
    C       = zeros(numr,numc);
    p       = randperm(numr*numc);
    L       = round(outlier_den*numr*numc);
    mgB     = max(abs(X0(:)));
    C(p(1:L))  = outlier_fac * mgB * rand(L,1);
    
    X = X + C;
    XzeroOutliers       = X;
    XzeroOutliers(C~=0) = 0;
    
    X_vec = X(:)';
    X_miss = X_vec.*ind_vec;
    X_miss(X_miss == 0) = [];
    
    X_matrix_sparse = sparse(I,J,X_miss,numr,numc);
    Indicator = sparse(I,J,1,numr,numc);
    Indicator_full = full(Indicator); 
    X_matrix_missing = full(X_matrix_sparse);
    X = X_matrix_missing;
    
 
    %% Algorithms

    OPTS = [];
    fprintf('+ NORST: ')
    t_start = tic;
    [~,PER]   = NORST_ST(X,V,Indicator_full);
    toc(t_start)
    rho_1 = rho_1 + PER.SEP;
   
  
end

rho_1 = rho_1/(r*nb_exp);



makerSize = 11;
numbMarkers = 500;
LineWidth = 2;

color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];

blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];


%% SEP
fig = figure;
k = numbMarkers;
d1 = semilogy(1:k:numc,rho_1(1:k:end),...
    'marker','d','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);
lgd = legend([d1],'NORST');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('SEP','interpreter','latex','FontSize',13,'FontName','Times New Roman');

set(fig, 'units', 'inches', 'position', [0.5 0.5 7 6]);
h=gca;
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:1000:5000,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
  'FontName','Times New Roman');
set(h,'FontSize', 22);
axis([0 numc 1e-6 1e2]); 
grid on;

