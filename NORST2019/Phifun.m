function y = Phifun(x,flag,idx,P)
    n = size(P,1);
    I = speye(n);
    if flag == 1    %A*x
        y = I(:,idx) * x - P *(P(idx,:)'*x);
    elseif flag == 2    %A'*x
        y = I(idx,:) * x - P(idx,:) * (P' * x);
    else
        error("Error: flag should be 1 or 2\n")
    end
end