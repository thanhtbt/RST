function er = sub_est_per(U_tr,U_es,type)

    if nargin < 3
        type = 'SEP';
    end
    [n,r] = size(U_tr);
    switch type
        case 'SEP'
       		 % Nguyen, Viet-Dung, Karim Abed-Meraim, Nguyen Linh-Trung, and Rodolphe Weber.
       		 % "Generalized minimum noise subspace for array processing."
        	 %  IEEE Transactions on Signal Processing 65, no. 14 (2017): 3789-3802.
            U_tr = orth(U_tr);
            U_es = orth(U_es);
            er   = abs(trace(U_es' * (eye(n)-U_tr*U_tr') * U_es)/trace(U_es'*(U_tr*U_tr')*U_es));
        case 'RE'
            ER = U_tr - U_es * (U_es'*U_tr); 
            er = norm(ER,'fro');
        case 'Angle'
            er = subspace(U_tr,U_es);  
            er = sin(er);   
        case 'SE' 
            % Narayanamurthy, Praneeth, and Namrata Vaswani. 
            % "Provable dynamic robust PCA or robust subspace tracking." 
            % IEEE Transactions on Information Theory 65.3 (2018): 1547-1577.
            [U_tr, ~] = qr(U_tr);
            [U_es, ~] = qr(U_es); 
            er = norm((eye(n) - U_tr(:, 1 : r)  * U_tr(:, 1 : r)') * U_es(:, 1 : r), 2);
        case 'EV'
            er =  abs(trace(U_tr' * U_es * U_es' * U_tr) /  trace(U_es * U_es')); 
    end

end


