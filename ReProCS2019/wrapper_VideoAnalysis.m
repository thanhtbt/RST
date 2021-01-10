% This wrapper calls the NORST function to perform the task of Background
% Recovery for real data (videos) in presence of occlusions and corruptions

clear;
clc;


addpath('YALL1_v1.4') % ell-1 minimization package
addpath('data') % contains of sample videos

%% Loading video data matrix
video = ["Curtain","SwitchLight","Lobby"];
PATH = '/home/vahidd/Git/NORST-rmc/data';
load([PATH,'/',char(video(1)),'.mat'])

L = I;  % video with foreground
Train = DataTrain;  % video with background only

%% Parameter Initialization
[n,m] = size(L);
r = 30;
t_max = m;
alpha = 60;

%% Generating missing entries' support
% (a)-Bernoulli Model:

% rho = 0.1; %denotes fraction of missing entries
% BernMat = rand(n, t_max);
% T = 1 .* (BernMat <= 1 - rho);
    
% (b)-Moving Model: (rectangular object moving along the width of the video)
T = ones(n,m);

height = imSize(1);
width = imSize(2);    

% b-a = rectangle's width
a = floor(height/2) - 2;
b = a + 4;

% frames with moving object
idx_frame = [width * 0 + 1 : width * floor(t_max/width)];

smin = 0;
smax = 1;
for j = idx_frame   
    for i = smin:smax
        T(height*i+ a : height*i + b ,j) = zeros(b-a+1,1);
    end
    smax = smax+1;
    if(smax - smin > width/4)
        smin = smin + 1;
    end

    if(smax >= width)
        smax = 1;
        smin = 0;
    end
end

M = L .* T; % corrupted version of the video frames
    
%% Calling NORST
fprintf('NORST\n')

% algorithm parameters
K = 3;
ev_thresh = 2e-3;
omega = 15;
tol = 1e-3; % tolerance in cgls and ncrpca functions

% mean subtraction
mu = mean(Train,2);
M_norst = M - mu;

t_norst = tic;
% initialization of true subspace
fprintf('Initialization...\t');
P_init = orth(ncrpca(Train, r, tol, 100));
fprintf('Subspace initialized\n');

[FG,BG] = NORST_video(M_norst, mu, T, P_init, ev_thresh, alpha, K, omega,tol);

t_NORST = toc(t_norst);                

%% Display the reconstructed video
DisplayVideo(L, T, M, BG, imSize,'recovered_BackGround_video.avi')
