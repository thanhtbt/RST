function [U_est,PER] = GRASTA(X,U_tr,OPTS)

[n, N] = size(X);
r      = size(U_tr,2);

%=======================GRASTA parameters==================================
OPTIONS.QUIET               = 1;     % suppress the debug information

OPTIONS.MAX_LEVEL           = 20;    % For multi-level step-size,
OPTIONS.MAX_MU              = 15;    % For multi-level step-size
OPTIONS.MIN_MU              = 1;     % For multi-level step-size

OPTIONS.DIM_M               = n;  % your data's ambient dimension
OPTIONS.RANK                = r; % give your estimated rank

OPTIONS.ITER_MIN            = 200;    % the min iteration allowed for ADMM at the beginning
OPTIONS.ITER_MAX            = 200;    % the max iteration allowed for ADMM
OPTIONS.rho                 = 2;      % ADMM penalty parameter for acclerated convergence
OPTIONS.TOL                 = 1e-8;   % ADMM convergence tolerance

OPTIONS.USE_MEX             = 0;     % If you do not have the mex-version of Alg 2
                                     % please set Use_mex = 0.
                                     
CONVERGE_LEVLE              = 20;    % If status.level >= CONVERGE_LEVLE, robust mc converges

OPTS  = struct(); % initial a empty struct for OPTS
Ut    = zeros(1); % U_hat will be initialized in GRASTA
status = struct();  % maintain GRASTA running status
status.init  = 0;   % will be set 1 once GRASTA start working
%==========================================================================

% performance assessment
PER.SEP     = zeros(1,N);
PER.Angle   = zeros(1,N);
PER.EV      = zeros(1,N);
PER.SE      = zeros(1,N);


%Processing
for k = 1:N
    x = X(:,k);
    idx = find(x); 
    % idx = (1:n)';% work with full data
    %====================GRASTA algorithm==================================
    [Ut, status, OPTS] = grasta_stream(x(idx), idx, Ut, status, OPTIONS, OPTS);
    if status.level >= CONVERGE_LEVLE,
        %         fprintf('Pass %d/%d, reach the convergence level - %d...\n',outIter, maxCycles,status.level);
        break;
    end
    

    PER.SEP(1,k)   = sub_est_per(U_tr,Ut,'SEP');
    PER.Angle(1,k) = sub_est_per(U_tr,Ut,'Angle');
    PER.SE(1,k)    = sub_est_per(U_tr,Ut,'SE');
    PER.EV(1,k)    = sub_est_per(U_tr,Ut,'EV');
   
end
U_est = Ut;
end
