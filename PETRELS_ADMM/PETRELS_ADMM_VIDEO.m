function [Wt,Rinv,bg_img,fg_img] = PETRELS_ADMM_VIDEO(x,idx,Wt,Rinv,OPTS)

%% Algorithm PETRELS-ADMM for subspace tracking in the presence of missing data and outliers
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
   

if isfield(OPTS,'lambda'),
    lambda = OPTS.lambda;
else
    lambda = 1;
end

[numr, r] = size(Wt);
x_Omega = x(idx,1);
W_Omega = Wt(idx,:);

OPTS = [];
s    = robustST_stream(x_Omega, W_Omega, OPTS);
indx = find(s==0);
x_re_free  = length(indx)/numr * x_Omega(indx);
epsilon    = 0;
if length(x_re_free) / numr  > epsilon
    W_Om_fre = W_Omega(indx,:);
    coeff = (W_Om_fre) \ x_re_free;
    p = W_Om_fre * coeff ;
    residual = x_re_free - p;
    eta = 1;
    for ii=1:length(indx)
        Tinv = Rinv(:,(idx(indx(ii))-1)*r+1:(idx(indx(ii)))*r);
        ss = Tinv*coeff;
        Tinv = (Tinv-ss*ss'*1/(lambda+coeff'*Tinv*coeff));
        Wt(idx(indx(ii)),:) =  W_Om_fre(ii,:) +  eta * lambda^(-1)*residual(ii)*coeff'*Tinv;
        Rinv(:,(idx(indx(ii))-1)*r+1:(idx(indx(ii)))*r) = Tinv;
    end
    Rinv = lambda^(-1)*Rinv;
else
    
end

W_Omega = Wt(idx,:);
w  =  Wt(idx(indx),:) \ x_re_free;  
bg_img = Wt * w; 
fg_img = x - bg_img;

noise_thres = min(x);
fg_img(abs(fg_img)<noise_thres) =0; 


end


