function [L,SEP,EER] = ROSETA_based_Subspace(YL,I,J,S, numr,numc, truerank, maxCycles,Params)

values = sparse(I,J,S,numr,numc);
Indicator = sparse(I,J,1,numr,numc);


%% ROSETA
lambda = 1/Params.HIGH_TEMP;

L_prev = orth(randn(numr,truerank));

W_ext = orth(YL);
for outiter = 1:maxCycles
    col_order = randperm(numc);
    
    for k=1:numc,
        idx = find(Indicator(:,col_order(k)));
        idxc = find(~Indicator(:,col_order(k)));
        Mask = ones(numr,1);
        Mask(idxc) = 0;
        v_Omega = zeros(numr,1);
        v_Omega(idx) = values(idx,col_order(k));
        
        
        
        [L, R, S, ParamsNew] = ROSETA(Mask, v_Omega, L_prev, truerank, Params,lambda);
        mu_adap(k) = ParamsNew.mu;
        Params = ParamsNew;
        L_prev = L;
     
        
        
        TS = L' * (eye(numr) - W_ext*W_ext') * L;
        MS = L' * (W_ext*W_ext') * L;
        SEP(outiter, k) = trace(TS)/trace(MS);
        EER = 0;
%         TS = L*L' - W_ext*W_ext';
%         TS = norm(TS,'fro');
%         MS = W_ext*W_ext';
%         MS = norm(MS,'fro');
%         EER(outiter,k) = TS/MS;
        

    end
end