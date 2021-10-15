
%% First we generate the data matrix and incomplete sample vector.
clc; close all;
clear all;


outlierFac = 0;
SAMPLING   = 0.9;
noiseFac = 1e-6;
nMonteCarlo = 2; % >=2

% Number of rows and columns
numr = 200;
numc = 1000;    
probSize = [numr,numc];
truerank = 6;
N = numr*numc;
M = round(SAMPLING * N);

% ROSETA
Params.const = 10; % parameter for the adaptive stepsize update
Params.val = 10; % decides how fast the stepsize (\mu) can change
Params.temp = 10; % initial hyper parameter
Params.HIGH_TEMP = 500;
Params.mu = max(Params.const./(1+Params.temp), 1); % initial \mu
Params.mufloor = 0.03;
mu_adap = zeros(numc,1); % save the computed adaptive \mu

% GRASTA
OPTIONS.QUIET               = 1;     % suppress the debug information

OPTIONS.MAX_LEVEL           = 20;    % For multi-level step-size,
OPTIONS.MAX_MU              = 15;    % For multi-level step-size
OPTIONS.MIN_MU              = 1;     % For multi-level step-size
OPTIONS.DIM_M               = numr;  % your data's ambient dimension
OPTIONS.RANK                = truerank; % give your estimated rank
OPTIONS.ITER_MIN            = 20;    % the min iteration allowed for ADMM at the beginning
OPTIONS.ITER_MAX            = 20;    % the max iteration allowed for ADMM
OPTIONS.rho                 = 2;   % ADMM penalty parameter for acclerated convergence
OPTIONS.TOL                 = 1e-8;   % ADMM convergence tolerance
OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.                                    
CONVERGE_LEVLE              = 20;    % If status.level >= CONVERGE_LEVLE, robust mc converges


if nMonteCarlo < 2
    nMonteCarlo = 2;
end
for iTest = 1 : nMonteCarlo,

fprintf('Number of Test: %d\n',iTest)
% The left and right factors which make up our true data matrix Y.
YL = randn(numr,truerank);
W_ext = orth(YL);
YR = randn(numc,truerank);

% Select a random set of M entries of Y.
idx = unique(ceil(N*rand(1,(10*M))));
idx = idx(randperm(length(idx)));

[I,J] = ind2sub([numr,numc],idx(1:M));
[J, inxs]=sort(J'); I=I(inxs)';

% Values of Y at the locations indexed by I and J.
S = sum(W_ext(I,:).*YR(J,:),2);
% S = sum(YL(I,:).*YR(J,:),2);
S_noiseFree = S;


noise = noiseFac*max(S)*randn(size(S));
S = S + noise;

% Add sparse outliers
outlier_magnitude = 0.01 * max(abs(S));
idx = randperm(length(S));
sparseIdx = idx(1:ceil(outlierFac*length(S)))';
Outlier_part =   randn(size(sparseIdx));  % randi([-1000, 1000],size(sparseIdx)); 

S(sparseIdx) = S(sparseIdx) +  outlier_magnitude*Outlier_part;
maxCycles  =1;

% ROSETA

disp('Update U(t) - ROSETA');
[U_roseta,sep,err] = ROSETA_based_Subspace(YL,I,J,S, numr,numc, truerank, maxCycles,Params);
SEP(iTest,:) = sep;
err_reg(iTest,:) = err;



end


SEP = sum(SEP)/nMonteCarlo;
% SEP1 = sum(SEP1)/nMonteCarlo;



%%  Check subspace 
figure(1);
k=1;
semilogy(1:k:numc,SEP(1:k:numc),'r-','MarkerSize',10,'linewidth',1.5); hold on;
legend('ROSETA','GRASTA');
grid on;
xlabel('data stream index');
ylabel('SEP');

%% Check outliers
% % % % % % % % % % % % % % % % % % % % % %
% Random select one column to show the recovery result.
% % % % % % % % % % % % % % % % % % % % % %

% p = randperm(numc); iSelect = p(1);
% S_free = YL*YR';
% Col_data = S_free(:,iSelect);
% 
% Col_max = max(abs(Col_data));
% idx = randperm(numr); sparseIdx = idx(1:ceil(outlierFac*numr))';
% Col_outlier = zeros(numr,1);
% Col_outlier(sparseIdx) = 1 * Col_max * randn(size(sparseIdx));
% Col_noise = noiseFac * Col_max * randn(numr,1);
% Col_data = Col_data + Col_outlier + Col_noise;
% 
% OPTS2.MAX_ITER      = 30;
% OPTS2.TOL           = 1e-8;


