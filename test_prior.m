rng default; % For reproducibility
X = [randn(100,2)*0.75+ones(100,2);
    randn(100,2)*0.5-ones(100,2)];
figure;
plot(X(:,1),X(:,2),'.');
title 'Randomly Generated Data';
[N,D] = size(X);
%%
K = 3;
clst_rej = [1,3];
clst_acc = 2;
theta_new = {1,3,K};
priors = ones(1,K)/K;
mu = X(randperm(N,K),:);
co_var_mat = X'*X;
S = 1;
for k = 1:K
   theta_new{S,1,k} = co_var_mat; %co variance matrix is stored here
   theta_new{S,2,k} = mu(k,:); % the mean for a cluster
   theta_new{S,3,k} = priors(1,k); % prior for a cluster
end
co_var_mats = {K};
for k = 1:K
   co_var_mats{k} = co_var_mat; %co variance matrix is stored here
   %theta_new{S,2,k} = mu(k,:); % the mean for a cluster
   %theta_new{S,3,k} = priors(1,k); % prior for a cluster
end
%q = zeros(K,N);
q = P_h_givn_x(X,k,priors,mu,co_var_mats)';
%q = repmat(ones(1,N)/K,K,1);
%disp(q);
%Q = [r,q(:,2:N)];
I_q_theta_thetaSNew( 1, K, q, theta_new , X, clst_rej, clst_acc )