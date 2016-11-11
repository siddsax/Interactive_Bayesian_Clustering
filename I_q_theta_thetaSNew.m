function [ sum_IQ ] = I_q_theta_thetaSNew( S, K, q, theta_old , X, clst_rej, clst_acc )
sum_IQ = 0;
[~,D] = size(X);
co_var_mat_s = {K};
mu_s = zeros(K,D);
priors_s = zeros(1,K);
for s = 1:S
   for k = 1:K
      co_var_mat_s{k} = theta_old{s,1,k};
      mu_s(k,:) = theta_old{s,2,k};
      priors_s(1,k) = theta_old{s,3,k};
   end 
   sum_IQ = sum_IQ + I_q_theta_thetaS( s, K, q,co_var_mat_s, priors_s, mu_s, X, clst_rej, clst_acc );
end
end
