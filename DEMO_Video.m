
clear; clc; 
run_path;

load('X_Lobby.mat');


%% Parameters
Rank                = 5;  % the estimated low-rank

disp(' start ... !')
t_start = tic;
[video_fg,video_bg,video_matrix] = Video_BGFG_Separation(X_video,Rank);
t_end = toc(t_start);
fprintf(' run time: %f (s) \n',t_end)


t = 363;
I_t = X_video(:,:,t);

BG  = video_bg(:,t); BG = reshape(BG,size(I_t));
FG  = video_fg(:,t); FG = reshape(FG,size(I_t));
fig = figure;
subplot(131); imagesc(I_t); axis off; colormap gray;
title('Original Frame')
subplot(132); imagesc(BG); axis off; colormap gray;
title('Background')
subplot(133); imagesc((FG)); axis off; colormap gray;
title('Foreground')
set(fig, 'units', 'inches', 'position', [0.2 0.5 13 6]);



t = 641;
I_t = X_video(:,:,t);

BG  = video_bg(:,t); BG = reshape(BG,size(I_t));
FG  = video_fg(:,t); FG = reshape(FG,size(I_t));
fig = figure;
subplot(131); imagesc(I_t); axis off; colormap gray;
title('Original Frame')
subplot(132); imagesc(BG); axis off; colormap gray;
title('Background')
subplot(133); imagesc((FG)); axis off; colormap gray;
title('Foreground')
set(fig, 'units', 'inches', 'position', [0.2 0.5 13 6]);

t = 1246;
I_t = X_video(:,:,t);

BG  = video_bg(:,t); BG = reshape(BG,size(I_t));
FG  = video_fg(:,t); FG = reshape(FG,size(I_t));
fig = figure;
subplot(131); imagesc(I_t); axis off; colormap gray;
title('Original Frame')
subplot(132); imagesc(BG); axis off; colormap gray;
title('Background')
subplot(133); imagesc((FG)); axis off; colormap gray;
title('Foreground')
set(fig, 'units', 'inches', 'position', [0.2 0.5 13 6]);

