function [fgs_cmp,bgs_cmp,video_matrix] = mybgfg_seperation(Video,Rank)

%%  PETRELS-ADMM for video background/foreground seperation
% Author      : Le Trung Thanh
% Email       : letrungthanhtbt@gmail.com // thanhletrung@vnu.edu.vn
% Address     : Vietnam National Unviersity, Hanoi
%               University of Engineering and Technoglogy
%               707 E3 Building, 144 Xuan Thuy Road, Hanoi City, Vietnam


Video_Length = size(Video,3);
rows         = size(Video,1);
colms        = size(Video,2);
n            = rows * colms;

Z_hat           = repmat(100*eye(Rank),1, n);
fgs_cmp         = zeros(n, Video_Length);
bgs_cmp         = zeros(n, Video_Length);
video_matrix    = zeros(n, Video_Length);


U_hat           = orth(randn(n,Rank));

for ii = 1 : Video_Length
    if mod(ii,20) == 0 
        fprintf('frame %d\n',ii)
    end
      
    I    = Video(:,:,ii);
    I    = I/max(max(I));
    I_Omega = I(:);
    idx = 1:length(I_Omega); % find(I_Omega);  % missing;
    % idxt = find(~I_Omega); % non-missing
    OPT = [];
    [U_hat,Z_hat,bg_img,fg_img] = PETRELS_ADMM_VIDEO(I_Omega,idx,U_hat,Z_hat,OPT);
    
    fgs_cmp(:,ii) = fg_img;
    bgs_cmp(:,ii) = bg_img; 
    
    
end
end

