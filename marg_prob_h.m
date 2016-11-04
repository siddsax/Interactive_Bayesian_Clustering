function [ P_h ] = marg_prob_h( X, co_var_mat_h, prior_h, mu_h )
%[N,~] = size(X); 
co_var_mat{1} = co_var_mat_h;
prior(1,1) = prior_h;
mu(1,:) = mu_h;
P_h = sum(P_h_givn_x(X, 1, prior, mu, co_var_mat));    
end


