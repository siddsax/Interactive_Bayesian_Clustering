% theta is having the means, co_var_mat, prior{pi_k in notes} for one iteration
% the inputs are 1) X is the data,2) max_iter is the number of iteratins for
% EM, 3) K is the number of clusters 3) theta_old is cell having the thetas
% of the previous steps 4) clst_rej are the clusters that were rejected in
% the last step 5) clst_acc are the clusters that were accepted in
% the last step 
function [ theta , theta_new ] = EM( X, max_iter, K, theta_old, clst_rej, clst_acc )
[N,~] = size(X);
mu = X( randsample(N,K),:);
co_var_mat = cell(K);
priors = ones(1,K)/K;
q = repmat(ones(1,N)/K,K,1);
maxIterCD = 500;
if nargin == 6
   [S,~,~]=size(theta_old);  %Calculate the number of previous iterations
end
co_var_mat_init = transpose(X)*X;
disp(det(co_var_mat_init));
for i = 1:K
    co_var_mat{i} = co_var_mat_init;
end
N_ks = zeros(1,K);
iter = 0;
log_lik = [];
disp('All Variable initialized');
while iter < max_iter
   disp(iter);
   for i = 1:K
   %The E step
       %Calculate P(h|theta, x_j) for all j 
       P_h_given_x = transpose(P_h_givn_x(X, K, priors, mu, co_var_mat));
       disp('geez');
       %implement stochastic co-ordinate descent using I_q_theta_thetaS and
       %KLDiv to get the 'q'(arbitrary probability) matrix
       if nargin == 6
          iterCD = 0;
          qTemp = zeros(K,N);
          while (iterCD < maxIterCD) && sum(sum(abs(q-qTemp))) < epsilon
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
       elseif nargin == 3 
           q = P_h_given_x;
       else
           disp('Wrong number of arguments in EM');
       end
       disp('E Step Done');
   %the M Step
       % Update steps as a normal GMM
       for j = 1:N
          N_ks(1,i) = N_ks(1,i) + q(i,j);%kXn
          mu(i,:) = mu(i,:) + q(i,j)*X(j,:);
          co_var_mat{i} = co_var_mat{i} + q(i,j)*(X(i,:) - mu(i,:))*(X(i,:) - mu(i,:)).';
       end
       disp('M Step Done');
       mu(i,:) = (1/N_ks(1,i))*mu(i,:);
       co_var_mat{i} = (1/N_ks(1,i))*co_var_mat{i};
       priors(1,i) = N_ks(1,i)/N;
   end
   P_h_given_x = P_h_givn_x(X, K, priors, mu, co_var_mat);
   log_lik_iter = sum(prior*transpose(P_h_given_x));
   % finding the log likelihood i.e $\sum_{i = 1}^{N}\sum_{j = 1}^{K}\pi_k*log(P(h,x|\theta)$
   log_lik = [log_lik,log_lik_iter];
   disp(log_lik);
   disp('Log Likelihood updated');
   if iter > 2
       delta = log_lik(1,-1) - log_lik(1,-2);
       if delta < epsilon
           disp('EM has converged');
           break;
       end
   end
   iter = iter + 1;
end
theta = {K,3};
% returning thr parameters of the model in theta which is a cell of size
% {K,3}
for k = 1:K
    theta{k,1} = co_var_mat{k}; %co variance matrix is stored here
    theta{k,2} = mu(k,:); % the mean for a cluster
    theta{k,3} = priors(1,k); % prior for a cluster
end
%append the now generated theta at the end of theta_old
if nargin == 6
    theta_old{S+1} = theta;
    theta_new = theta_old;
else
    theta_new{1} = theta;
end    
end
