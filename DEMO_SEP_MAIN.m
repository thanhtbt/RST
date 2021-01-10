%% DEMO: Robust Subspace Tracking with Missing data and Outliers
% Author      : Le Trung Thanh
% Email       : letrungthanhtbt@gmail.com // thanhletrung@vnu.edu.vn
% Address     : Vietnam National Unviersity, Hanoi
%               University of Engineering and Technoglogy
%               707 E3 Building, 144 Xuan Thuy Road, Hanoi City, Vietnam

% Reference   : [1] L.T. Thanh, V-D. Nguyen, N.L. Trung, K. Abed-Meraim
%                   "Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee". 
%                   IEEE Trans. Signal Process. 2021 (accepted).
%               [2] L.T. Thanh, V-D. Nguyen, N. L. Trung and K. Abed-Meraim. 
%                   "Robust Subspace Tracking with Missing Data and Outliers via ADMM". 
%                   EUSIPCO, 2019.
%%
clear; clc;
run_path;


nb_exp            = 3;          % Number of independent runs
outlier_intensity = 10;         % Outlier Intensity (Magnitude)
outlier_density   = 0.2;        % Outlier density 
SAMPLING          = 0.8;        % Observations  
SNR               = 20;         % dB
noise_level       = 1 * 10^(-SNR/20); 

numr = 50;    % Number of rows
numc = 2001;  % Number of columns
r = 5;        % Rank

% Evaluation Metrics
Angle_ROSETA  = zeros(1,numc);      SEP_ROSETA  = zeros(1,numc);
Angle_GRASTA  = zeros(1,numc);      SEP_GRASTA  = zeros(1,numc);
Angle_ReProCS = zeros(1,numc);      SEP_ReProCS = zeros(1,numc);
Angle_NORST   = zeros(1,numc);      SEP_NORST   = zeros(1,numc);
Angle_OURS    = zeros(1,numc);      SEP_OURS    = zeros(1,numc);

SE_ROSETA  = zeros(1,numc);
SE_GRASTA  = zeros(1,numc);
SE_ReProCS = zeros(1,numc);
SE_NORST   = zeros(1,numc);
SE_OURS    = zeros(1,numc);


% ============================== PROGRAM  ==============================

for ii_exp = 1: nb_exp
    
    fprintf('\nExperience %d of %d ===========================================\n', ii_exp,nb_exp)
   
    %% Data Generation
    S       = randn(r,numc);
    A       = randn(numr,r);
    X0      = A*S;
    V       = orth(X0);  V = V(:,1:r);
    normX   = norm(X0,'fro');
    
    % Gaussian Noise
    N      = randn(numr, numc);
    normN  = norm(N,'fro');
    X      = X0 + noise_level * normX/normN * N;
    X_free = X;
    
    % Outliers
    C   = zeros(numr,numc);
    p   = randperm(numr*numc);
    L   = round(outlier_density*numr*numc);
    mgB = max(abs(X0(:)));
    C(p(1:L))  = outlier_intensity * mgB * rand(L,1);
    
    % Missing + Outliers 
    
    ind      = rand(numr,numc);
    Mask = 1 .* (ind <= SAMPLING);
    X         = X + C;
    X         = X .* Mask;  
    
    
    %% State-of-the-art RST Algorithms
    % ======================== ROSETA ====================================%
    t_start = tic;
    [~,PER_ROSETA] = ROSETA_ST(X,V);
    t_end = toc(t_start); 
    SEP_ROSETA     = SEP_ROSETA   + PER_ROSETA.SEP;
    Angle_ROSETA   = Angle_ROSETA + PER_ROSETA.Angle;
    SE_ROSETA      = SE_ROSETA    + PER_ROSETA.SE;
    fprintf('+ ROSETA-ST: %f (s) \n',t_end);
    
    % ======================== GRASTA ====================================%
   
    t_start = tic;
    [~,PER_GRASTA] = GRASTA(X,V);
    t_end = toc(t_start);
    SEP_GRASTA   = SEP_GRASTA   + PER_GRASTA.SEP;
    Angle_GRASTA = Angle_GRASTA + PER_GRASTA.Angle;
    SE_GRASTA    = SE_GRASTA    + PER_GRASTA.SE;
    fprintf('+ GRASTA: %f (s)\n',t_end);
    
    
    % ======================== ReProCS ====================================%
  
    OPTS = [];
    t_start = tic;
    [~,PER_ReProCS] = ReProCS_ST(X,V,OPTS);
    t_end = toc(t_start);
    SEP_ReProCS   = SEP_ReProCS   + PER_ReProCS.SEP;
    Angle_ReProCS = Angle_ReProCS + PER_ReProCS.Angle;
    SE_ReProCS    = SE_ReProCS    + PER_ReProCS.SE;
    fprintf('+ ReProCS: %f (s)\n',t_end)
    
    
    t_start = tic;
    [~,PER_NORST]   = NORST_ST(X,V,ind);
    t_end = toc(t_start);
    SEP_NORST   = SEP_NORST    + PER_NORST.SEP;
    Angle_NORST = Angle_NORST  + PER_NORST.Angle;
    SE_NORST    = SE_NORST     + PER_NORST.SE;
    fprintf('+ NORST: %f (s)\n',t_end)
    % ======================== ADMM ====================================%
   

    OPTS = [];
    OPTS.method  = 'ADMM';
    t_start      = tic;  
    [~,PER_ADMM] = PETRELS_ADMM(X,V,OPTS);
    t_end = toc(t_start);
    SEP_OURS   = SEP_OURS  + PER_ADMM.SEP;
    Angle_OURS = Angle_OURS+ PER_ADMM.Angle;
    SE_OURS    = SE_OURS   + PER_ADMM.SE;
    fprintf('+ PETRELS-ADMM: %f (s) \n',t_end);

    
end

Angle_ROSETA  = Angle_ROSETA/nb_exp;     SEP_ROSETA  = SEP_ROSETA/(nb_exp);
Angle_GRASTA  = Angle_GRASTA/nb_exp;     SEP_GRASTA  = SEP_GRASTA/(nb_exp);
Angle_ReProCS = Angle_ReProCS/nb_exp;    SEP_ReProCS = SEP_ReProCS/(nb_exp);
Angle_NORST   = Angle_NORST/nb_exp;      SEP_NORST   = SEP_NORST/(nb_exp);
Angle_OURS    = Angle_OURS/nb_exp;       SEP_OURS    = SEP_OURS/(nb_exp);

SE_ROSETA = SE_ROSETA/(nb_exp);
SE_GRASTA = SE_GRASTA/(nb_exp);
SE_ReProCS = SE_ReProCS/(nb_exp);
SE_NORST = SE_NORST/(nb_exp);
SE_OURS= SE_OURS/(nb_exp);


%% plot

makerSize = 11;
numbMarkers = 500;
LineWidth = 2;
k = round(numc/10);

color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
magenta_0 = [1 0 1];
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
subplot(131);

d1 = semilogy(1:k:numc,SEP_ROSETA(1:k:end),...
    'marker','d','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);
hold on;
d2 = semilogy(1:k:numc,SEP_GRASTA(1:k:end),...
    'marker','^','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
hold on;
d3 = semilogy(1:k:numc,SEP_ReProCS(1:k:end),...
    'marker','o','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
hold on;
d4 = semilogy(1:k:numc,SEP_NORST(1:k:end),...
    'marker','*','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',brow_n,'LineWidth',LineWidth);
hold on;
d10 = semilogy(1:k:numc,SEP_OURS(1:k:end),...
    'marker','s','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
hold on;

lgd = legend([d1,d2,d3,d4,d10],'ROSETA','GRASTA','ReProCS','NORST','Proposed');
lgd.FontSize = 14;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('SEP($\mathbf{U_{tr}},\mathbf{U_{es}}$)','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h1=gca;
set(h1,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h1,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
   'FontName','Times New Roman');
set(h1,'FontSize', 16);
axis([0 numc 1e-6 1e2]); 
grid on;

%% SE
subplot(132); 
d1 = semilogy(1:k:numc,SE_ROSETA(1:k:end),...
    'marker','d','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);
hold on;
d2 = semilogy(1:k:numc,SE_GRASTA(1:k:end),...
    'marker','^','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
hold on;
d3 = semilogy(1:k:numc,SE_ReProCS(1:k:end),...
    'marker','o','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
hold on;
d4 = semilogy(1:k:numc,SE_NORST(1:k:end),...
    'marker','*','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',brow_n,'LineWidth',LineWidth);
hold on;
d10 = semilogy(1:k:numc,SE_OURS(1:k:end),...
    'marker','s','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
hold on;

xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('SE($\mathbf{U_{tr}},\mathbf{U_{es}}$)','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h2 = gca;
set(h2,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h2,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
  'FontName','Times New Roman');
set(h2,'FontSize', 16);
axis([0 numc 1e-3 1e0]); 
grid on;


%% Angle
subplot(133);
d1 = semilogy(1:k:numc,Angle_ROSETA(1:k:end),...
    'marker','d','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);
hold on;
d2 = semilogy(1:k:numc,Angle_GRASTA(1:k:end),...
    'marker','^','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
hold on;
d3 = semilogy(1:k:numc,Angle_ReProCS(1:k:end),...
    'marker','o','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
hold on;
d4 = semilogy(1:k:numc,Angle_NORST(1:k:end),...
    'marker','*','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',brow_n,'LineWidth',LineWidth);
hold on;
d10 = semilogy(1:k:numc,Angle_OURS(1:k:end),...
    'marker','s','markersize',makerSize,'markerfacecolor','w',...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
hold on;

xlabel('Time Index','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\sin\big(\theta(\mathbf{U_{tr}},\mathbf{U_{es}})\big)$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
h3=gca;
set(h3,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h3,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
  'FontName','Times New Roman');
set(h3,'FontSize', 16);
axis([0 numc 1e-3 1e0]); 
grid on;

set(fig, 'units', 'inches', 'position', [0.2 0.5 13 6]);


