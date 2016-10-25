function [ req_probab ] = P_h_hs( q, co_var_mat_hs, prior_hs, X , mu_hs) 
[N,~] = size(X); 
co_var_mat{1} = co_var_mat_hs;
priors(1,1) = prior_hs;
mu(1,:) = mu_hs;
req_probab = q*P_h_givn_x(X, 1, priors, mu, co_var_mat); %nXk   
req_probab = req_probab/N;
end

