function [ theta_new,BIC,P_h_given_x,mu ] = gmm( X,K)

[N,D] = size(X);

iter = 1;
max_iter = 20; % Max iterations for EM

CLL = zeros(max_iter,1);
ILL = zeros(max_iter,1);

S = 1; %initializing the value for the meta iterations

co_var_mat = cell(K);   % initializing the  co variance matrix for gaussians
co_var_mat_init = transpose(X)*X;
for i = 1:K
    co_var_mat{i} = co_var_mat_init;
end

% initialize priors uniformly
priors = ones(1,K)/K; 

% we initialize the means totally randomly
mu = X( randsample(N,K),:); 

log_lik = zeros(1,max_iter);
while iter < max_iter
   disp(iter); 
   %The E step
   P_h_given_x = P_h_givn_x(X, K, priors, mu, co_var_mat); %NXK
   q = P_h_given_x; %NXK
   %disp('E Step Done');
   %pause;
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
BIC =  -2*log_lik(1,iter) + (K*D + D*D*K)*log(N);
end    

