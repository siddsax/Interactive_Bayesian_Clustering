% theta is having the means, co_var_mat, prior{pi_k in notes} for one iteration
% the inputs are 1) X is the data,2) max_iter is the number of iteratins for
% EM, 3) K is the number of clusters 3) theta_old is cell having the thetas
% of the previous steps 4) clst_rej are the clusters that were rejected in
% the last step 5) clst_acc are the clusters that were accepted in
% the last step 
function [ theta_new ] = EM( X, max_iter, K,epsil, S, theta_old, clst_rej, clst_acc )
[N,D] = size(X);

iter = 1;

CLL = zeros(max_iter,1);
ILL = zeros(max_iter,1);
BIC = zeros(max_iter,1);
co_var_mat = cell(K);   % initializing the  co variance matrix for gaussians
co_var_mat_init = transpose(X)*X;
for i = 1:K
    co_var_mat{i} = co_var_mat_init;
end

% initialize priors uniformly
priors = ones(1,K)/K; 

% we initialize the means totally randomly
mu = X( randsample(N,K),:);

q = repmat(ones(1,N)/K,K,1);
maxIterCD = 500;
beta = 1;
alpha = 100;
if nargin == 8
   [S,~,~]=size(theta_old);  %Calculate the number of previous iterations
end

log_lik = zeros(1,max_iter);
while iter < max_iter
   %disp(iter);
   %The E step
   %Calculate P(h|theta, x_j) for all j 
   P_h_given_x = P_h_givn_x(X, K, priors, mu, co_var_mat); %NXK
   %implement stochastic co-ordinate descent using I_q_theta_thetaS and
   %KLDiv to get the 'q'(arbitrary probability) matrix
   if nargin == 8
%      iterCD = 0;
      fun = @(q)beta*I_q_theta_thetaSNew( S, K, q, theta_old , X, clst_rej, clst_acc ) + alpha * KLDivNew(P_h_given_x,q);
      q = fminunc(fun,q);

%       while (iterCD < maxIterCD)
%           qTemp = q;
%           r=q(:,1);
%           fun = @(r)beta*I_q_theta_thetaSNew( S, K, [r,q(:,2:N)], theta_old , X, clst_rej, clst_acc ) + alpha * KLDivNew(P_h_given_x,[r,q(:,2:N)],K,N);
%           r = fminunc(fun,r);
%           q(:,1)=r;
%           for j = 2:N-1
%              r = q(:,j);
%              fun = @(r)beta*I_q_theta_thetaSNew( S, K, [q(:,1:j-1),r,q(:,j+1:N)], theta_old , X, clst_rej, clst_acc ) + alpha * KLDivNew(P_h_given_x,[q(:,1:j-1),r,q(:,j+1:N)],K,N);
%              r = fminunc(fun,r);
%              q(:,j)=r;
%           end
%           r = q(:,N);
%           fun = @(r)beta*I_q_theta_thetaSNew( S, K,  [q(:,1:N-1),r], theta_old , X, clst_rej, clst_acc ) + alpha * KLDiv(P_h_given_x,[q(:,1:N-1),r],K,N);
%           r = fminunc(fun,r);
%           q(:,N)=r;
%           iterCD = iterCD + 1;
%       end
   elseif nargin == 5
       q = P_h_given_x; %NXK
   else
       disp('Wrong number of arguments in EM');
   end
   %the M Step
   % Update steps as a normal GMM
   N_ks = sum(q,1);%1XK
   priors = N_ks/N; %1XK
   mu = P_h_given_x'*X; %KXD
   for k = 1:K mu(k,:) = mu(k,:)/N_ks(1,k); end %KXD
   for i = 1:K
       co_var_mat{i} = zeros(D,D);
       for j = 1:N
          co_var_mat{i} = co_var_mat{i} + q(j,i)*(X(j,:) - mu(i,:)).'*(X(j,:) - mu(i,:));
       end
       co_var_mat{i} = (1/N_ks(1,i))*co_var_mat{i};
   end
   %disp('M Step Done');
   % compute *complete* and *incomplete* log likelihoods
   cll = 0;
   ill = 0;
   deter = zeros(1,K);
   for k = 1:K
       deter(1,k) = det(co_var_mat{k});
   end    
   for n=1:N
     ill_tmp = zeros(1, K);
     for k=1:K
         ill_tmp(k) = ill_tmp(k) + log(priors(k)) - (0.5)*(X(n,:)-mu(k,:))*inv(co_var_mat{k})*(X(n,:)-mu(k,:))' - 0.5*D*log(2*pi) - 0.5*log(deter(1,k));
         cll = cll + P_h_given_x(n,k)*ill_tmp(k);
     end
     ill = ill + logsumexp(ill_tmp,2);
   end
   CLL(iter) = cll;
   ILL(iter) = ill;
   BIC(iter) =  -2*ill + (K*D + D*D*K)*log(N); 
   if iter > 5
      if BIC(iter-1) - BIC(iter) < .0001
          disp('EM has converged');
          break;
      end
   end
   iter = iter + 1;
   disp(BIC(iter-1));
end
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
