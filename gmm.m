function [ theta_new ] = gmm( X,K)
[N,D] = size(X);
S = 1;
mu = X( randsample(N,K),:);
co_var_mat = cell(K);
priors = ones(1,K)/K;
q = repmat(ones(1,N)/K,K,1);
co_var_mat_init = transpose(X)*X;
for i = 1:K
    co_var_mat{i} = co_var_mat_init;
end
iter = 1;
% disp('All Variable initialized');
% pause;
max_iter = 1000;
log_lik = zeros(1,max_iter);
while iter < max_iter
   %The E step
   P_h_given_x = transpose(P_h_givn_x(X, K, priors, mu, co_var_mat));
   q = P_h_given_x;
   disp('E Step Done');
   %pause;
   %the M Step
   N_ks = zeros(1,K);
   mu = zeros(K,D);
   for i = 1:K
       % Update steps as a normal GMM
       % Update the N_ks
       N_ks(1,i) = sum(q(i,:));%kXN
       for j = 1:N 
          mu(i,:) = mu(i,:) + q(i,j)*X(j,:);
       end
       mu(i,:) = (1/N_ks(1,i))*mu(i,:);
       co_var_mat{i} = zeros(D,D);
       for j = 1:N
          co_var_mat{i} = co_var_mat{i} + q(i,j)*(X(j,:) - mu(i,:)).'*(X(j,:) - mu(i,:));
       end
       co_var_mat{i} = (1/N_ks(1,i))*co_var_mat{i};
       priors(1,i) = N_ks(1,i)/N;
   end
   disp('M Step Done');
   %pause;
   P_h_given_x = P_h_givn_x(X, K, priors, mu, co_var_mat);
   log_lik_iter = sum(priors*transpose(P_h_given_x));
   % finding the log likelihood i.e $\sum_{i = 1}^{N}\sum_{j = 1}^{K}\pi_k*log(P(h,x|\theta)$
   log_lik(1,iter) = log_lik_iter;
   disp('Log Likelihood updated');
   disp(log_lik(1,iter));
   %pause;
   if iter > 3
       delta = log_lik(1,iter) - log_lik(1,iter-2);
       %disp(delta);
       if delta < 0
           disp('EM has converged');
           pause;
           break;
       end
   end
   iter = iter + 1;
end
theta_new = {S,K,3};    
% returning thr parameters of the model in theta which is a cell of size
% {K,3}
%append the now generated theta at the end of theta_old
for k = 1:K
    theta_new{S,1,k} = co_var_mat{k}; %co variance matrix is stored here
    theta_new{S,2,k} = mu(k,:); % the mean for a cluster
    theta_new{S,3,k} = priors(1,k); % prior for a cluster
end
end    

