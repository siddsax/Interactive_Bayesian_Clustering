function [ req_probab ] = P_h_hs( q, co_var_mat, priors, X, mu, K) 
[N,~] = size(X); 
% q is KXN
% P_h_givn_x is a NXK
req_probab = q*P_h_givn_x(X, K, priors, mu, co_var_mat); %NXK   
req_probab = req_probab/N;
end

