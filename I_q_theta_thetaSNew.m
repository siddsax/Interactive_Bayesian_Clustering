function [ sum_IQ ] = I_q_theta_thetaSNew( S, K, q, theta_old , X, clst_rej, clst_acc )
sum_IQ = 0;
co_var_mat_s = {K};
mu_s = {K};
priors_s = {K};
for s = 1:S
   for k = 1:K
      co_var_mat_s{k} = theta_old{s,k,1};
      disp(class(co_var_mat_s{k}));
      disp('jeez');
      mu_s{k} = theta_old{s,k,2};
      priors_s{k} = theta_old{s,k,3};
   end 
   sum_IQ = sum_IQ + I_q_theta_thetaS( s, K, q,co_var_mat_s, priors_s, mu_s, X, clst_rej, clst_acc );
end
end
