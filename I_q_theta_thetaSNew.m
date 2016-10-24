function [ sum_IQ ] = I_q_theta_thetaSNew( S, K, q, theta_old , X, clst_rej, clst_acc )
sum_IQ = 0;
  for s = 1:S
      for k = 1:K
       sum_IQ = sum_IQ + I_q_theta_thetaS( s, K, q,theta_old{s,k,1}, theta_old{s,k,3}, theta_old{s,k,2}, X, clst_rej, clst_acc );
      end
  end
end
