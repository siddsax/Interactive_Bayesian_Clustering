function [ req_probab ] = P_h_hs( q, co_var_mat_hs, prior_hs, X , mu_hs) 
req_probab = 0;
[N,~] = size(X); 
for i = 1:N
    req_probab = req_probab + q(1,i)*P_h_givn_x(i , X, co_var_mat_hs, prior_hs, mu_hs);
end    
req_probab = req_probab/N;
end

