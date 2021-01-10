function [U_est,PER,U_Var] = PETRELS_ADMM(X,U_tr,OPTS)

%% Algorithm PETRELS-ADMM for subspace tracking in the presence of missing data and outliers
% Author      : Le Trung Thanh
% Email       : letrungthanhtbt@gmail.com // thanhletrung@vnu.edu.vn
% Address     : Vietnam National Unviersity, Hanoi
%               University of Engineering and Technoglogy
%               707 E3 Building, 144 Xuan Thuy Road, Hanoi City, Vietnam

% Input:
%               X     : Measurement Matrix (contains outliers + missing
%                       entries, e.g., X(i,j) = 0); 
%               U_tr  : True subspace 
% Output:  
%               U_est : Estimated subspace 
%               PER   : Performance metrics (SEP,SE,Angle,EV)
%                       

% Reference   : [1] L.T. Thanh, V-D. Nguyen, N.L. Trung, K. Abed-Meraim
%                   "Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee". 
%                   IEEE Trans. Signal Process. 2021 (accepted).
%               [2] L.T. Thanh, V-D. Nguyen, N. L. Trung and K. Abed-Meraim. 
%                   "Robust Subspace Tracking with Missing Data and Outliers via ADMM". 
%                   EUSIPCO, 2019.

%% Initializations

if isfield(OPTS,'lambda'), % Forgetting factor
     lambda = OPTS.lambda;
else lambda = 0.9999;
end
if isfield(OPTS,'outlier_den'), % Outlier density 
    outlier_den = OPTS.outlier_den;
else
    outlier_den = 0.0;
end
if isfield(OPTS,'gamma_thres'), 
    gamma_thres = OPTS.gamma_thres;
else
    gamma_thres = sin(0); 
end
 
[numr,numc] = size(X);
r           = size(U_tr,2);

PER.EV    = zeros(1,numc);
PER.SEP   = zeros(1,numc);
PER.SE    = zeros(1,numc);
PER.Angle = zeros(1,numc);
U_Var     = zeros(1,numc);

%% Initialization
Ut     = randn(numr,r);
Rinv   = repmat(100*eye(r),1,numr);
Ut_Old = Ut; 


%% Main Algorithm 
for k = 1:numc, % 
    x = X(:,k);
    
    idx_data_free = find(x);
    miss_den      = length(idx_data_free)/numr;
    x_Omega       = x(idx_data_free);
    U_Omega       = Ut(idx_data_free,:);
   
    [s,~,~]          = robustST_stream(x_Omega, U_Omega,OPTS);
    idx_outlier_free = find(s==0);
    x_re_free        = length(idx_outlier_free)/numr * x_Omega(idx_outlier_free);
    epsilon          = 1 - miss_den - outlier_den;
    
    if (length(x_re_free) / numr)  >= epsilon
        U_Omega_fre = U_Omega(idx_outlier_free,:);
        coeff    = U_Omega_fre \ x_re_free; 
        p        = U_Omega_fre * coeff ;
        residual = x_re_free - p;
        theta    = (norm(residual) * norm(p));
        gamma    =  theta / (sqrt(theta^2 + 1)) ; 
        if gamma > gamma_thres  eta_t = 1;  else eta_t = 1/gamma;  end
        lambda_t =  lambda * eta_t;
        mu_t     = 0;
        for ii=1:length(idx_outlier_free)                
            Tinv = Rinv(:,(idx_data_free(idx_outlier_free(ii))-1)*r+1 ...
                                 :(idx_data_free(idx_outlier_free(ii)))*r);
            ss   = (Tinv + mu_t * eye(r))*coeff;
            Tinv = (Tinv-ss*ss'*1/(lambda_t + coeff'*Tinv*coeff));
            
            Ut(idx_data_free(idx_outlier_free(ii)),:) =  U_Omega_fre(ii,:) ...
                    + lambda_t^(-1) * miss_den * residual(ii)*coeff'*Tinv;
                
            Rinv(:,(idx_data_free(idx_outlier_free(ii))-1)*r+1 ... 
                           :(idx_data_free(idx_outlier_free(ii)))*r) = Tinv;
        end
        Rinv = lambda_t^(-1) * Rinv;
    else
    end
%%  Performance Estimation 
    PER.SEP(1,k)   = sub_est_per(U_tr,Ut,'SEP');
    PER.Angle(1,k) = sub_est_per(U_tr,Ut,'Angle');
    PER.SE(1,k)    = sub_est_per(U_tr,Ut,'SE');
    PER.EV(1,k)    = sub_est_per(U_tr,Ut,'EV');
    U_Var(1,k)     = norm(Ut - Ut_Old,'fro')^2;
    Ut_Old = Ut;
end
%% Outputs 
    U_est = Ut; 

end


