% theta is having the means, co_var_mat, prior{pi_k in notes} for one iteration
% the inputs are 1) X is the data,2) max_iter is the number of iteratins for
% EM, 3) K is the number of clusters 3) theta_old is cell having the thetas
% of the previous steps 4) clst_rej are the clusters that were rejected in
% the last step 5) clst_acc are the clusters that were accepted in
% the last step 
function [ theta_new ] = EM( X, max_iter, K,epsil, S, theta_old, clst_rej, clst_acc )
[N,D] = size(X);
mu = X( randsample(N,K),:);
co_var_mat = cell(K);
priors = ones(1,K)/K;
q = repmat(ones(1,N)/K,K,1);
maxIterCD = 500;
if nargin == 8
   [S,~,~]=size(theta_old);  %Calculate the number of previous iterations
end
co_var_mat_init = transpose(X)*X;
for i = 1:K
    co_var_mat{i} = co_var_mat_init;
end
iter = 1;
% disp('All Variable initialized');
% pause;
log_lik = zeros(1,max_iter);
while iter < max_iter
   disp(iter);
   %pause;
   %The E step
       %Calculate P(h|theta, x_j) for all j 
       P_h_given_x = transpose(P_h_givn_x(X, K, priors, mu, co_var_mat));
       %implement stochastic co-ordinate descent using I_q_theta_thetaS and
       %KLDiv to get the 'q'(arbitrary probability) matrix
   if nargin == 8
      iterCD = 0;
      qTemp = zeros(K,N);
      while (iterCD < maxIterCD) 
          qTemp = q;
          fun = @(x)beta*I_q_theta_thetaSNew( S, K, [q(:,1),q(:,2:N)], theta_old , X, clst_rej, clst_acc ) + alpha * KLDivNew(P_h_given_x,[q(:,1),q(:,2:N)]);
          q(:,1) = fminsearch(fun,q(:,1));
          for j = 2:N-1
             fun = @(x)beta*I_q_theta_thetaSNew( S, K, [q(:,1:j-1),q(:,j),q(:,j+1:N)], theta_old , X, clst_rej, clst_acc ) + alpha * KLDivNew(P_h_given_x,[q(:,1:j-1),q(:,j),q(:,j+1:N)]);
             q(:,j) = fminsearch(fun,q(:,j));
          end
          fun = @(x)beta*I_q_theta_thetaSNew( S, K,  [q(:,1:N-1),q(:,N)], theta_old , X, clst_rej, clst_acc ) + alpha * KLDiv(P_h_given_x,[q(:,1:N-1),q(:,N)]);
          q(:,N) = fminsearch(fun,q(:,N));
          iterCD = iterCD + 1;
      end
   elseif nargin == 5
       q = P_h_given_x;
   else
       disp('Wrong number of arguments in EM');
   end
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
   if isnan(log_lik_iter)
       disp(priors);
       pause;
       disp(P_h_given_x);
       pause;
   end
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
disp(log_lik);
theta_new = {S,K,3};    
% returning thr parameters of the model in theta which is a cell of size
% {K,3}
%append the now generated theta at the end of theta_old
if nargin == 8
    for s = 1:S-1
       theta_new{s} = theta_old{s};
    end
    for k = 1:K
       theta_new{S,1,k} = co_var_mat{k}; %co variance matrix is stored here
       theta_new{S,2,k} = mu(k,:); % the mean for a cluster
       theta_new{S,3,k} = priors(1,k); % prior for a cluster
    end
else
    for k = 1:K
       theta_new{S,1,k} = co_var_mat{k}; %co variance matrix is stored here
       theta_new{S,2,k} = mu(k,:); % the mean for a cluster
       theta_new{S,3,k} = priors(1,k); % prior for a cluster
    end
end    
end
