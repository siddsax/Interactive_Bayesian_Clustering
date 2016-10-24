function [ P_h ] = marg_prob_h( X, co_var_mat_h, prior_h, mu_h )
[N,~] = size(X); 
P_h = 0;
for i = 1:N
    P_h = P_h  + P_h_givn_x(i, X, co_var_mat_h, prior_h, mu_h);
end



