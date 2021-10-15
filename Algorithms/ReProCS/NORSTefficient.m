function [L_hat] = NORSTefficient(P_init,M, T)
[n, t_max] = size(M);
L_hat = zeros(n,t_max);
P_hat = P_init;
[~,r] = size(P_hat);
I = speye(n);
%% Main Algorithm Loop
for ii = 1 : t_max

    T_miss = find(T(:, ii) == 1);
    y_t = M(T_miss,ii);
        
    tol = 1e-16;
    phi_t = I(:,T_miss) * I(:,T_miss)' * P_hat;
    a_hat = cgls(phi_t(T_miss,:), y_t, 0, tol, 20);    
    L_hat(:, ii) = P_hat * a_hat;
end

end