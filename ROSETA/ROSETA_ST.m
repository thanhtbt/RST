function [U_es,PER] = ROSETA(X,U_tr)


[numr, numc] = size(X);
r      = min(size(U_tr));

% performance assessment
PER.SEP   =   zeros(1,numc);
PER.Angle =   zeros(1,numc);
PER.SE    =   zeros(1,numc);
PER.EV    =    zeros(1,numc);

% Parametors for ROSEV
Params.const = 10; % parameter for the adaptive stepsize update
Params.val   =   10; % decides how fast the stepsize (\mu) can change
Params.temp  =  10; % initial hyper parameter
Params.HIGH_TEMP = 500;
Params.mu = max(Params.const./(1+Params.temp), 1); % initial \mu
Params.mufloor = 0.03;
mu_adap = zeros(numc,1); % save the computed adaptive \mu
lambda = 1/Params.HIGH_TEMP;
alpha = 0.1;

%% ROSETA

Ut = orth(randn(numr,r));

for k = 1:numc
    x    = X(:,k);
    idx  = find(x);
    idxc = find(~x);
    Mask = ones(numr,1);
    Mask(idxc) = 0;
    
    [Ut,~,~, ParamsNew] = ROSETA(Mask, x, Ut, r,Params,lambda);
    mu_adap(k) = ParamsNew.mu;
    Params     = ParamsNew;
    
    PER.SEP(k)     =  sub_est_per(U_tr,Ut,'SEP');
    PER.Angle(k)   =  sub_est_per(U_tr,Ut,'Angle');
    PER.SE(k)      =  sub_est_per(U_tr,Ut,'SE');
    PER.EV(k)      =  sub_est_per(U_tr,Ut,'EV');
    
end
U_es = Ut;
end