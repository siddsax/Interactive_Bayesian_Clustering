% theta is having the means, co_var_mat, prior{pi_k in notes} for one iteration
% the inputs are 1) X is the data,2) max_iter is the number of iteratins for
% EM, 3) K is the number of clusters 3) theta_old is cell having the thetas
% of the previous steps 4) clst_rej are the clusters that were rejected in
% the last step 5) clst_acc are the clusters that were accepted in
% the last step
function [ theta_new,score ] = EM( X, max_iter, K,~, theta_old, clst_rej, clst_acc )
[N,D] = size(X);
iter = 1;
frac_to_iter = .05;
eps = 10^-19;
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
q = repmat(ones(N,1)/K,1,K);
beta = 1;
alpha = 10^(-19);
if nargin == 7
   [S,~,~]=size(theta_old);  %Calculate the number of previous iterations
else 
    S = 0;
end

while iter < max_iter
   disp(iter);
   %The E step
   %Calculate P(h|theta, x_j) for all j
   P_h_given_x = P_h_givn_x(X, K, priors, mu, co_var_mat); %NXK
   %implement stochastic co-ordinate descent using I_q_theta_thetaS and
   %KLDiv to get the 'q'(arbitrary probability) matrix
   if nargin == 7
%% Coordinate Descent %%
    KL_avg=0;
    IQ_avg=0;
    fval_avg=0;
    co_=0;
    prev_cd=0;
    new_cd=0;
    q = P_h_given_x; %NXK
    q = q';
    epsilon_val = -1;
    while(co_< 20 && epsilon_val ~= 0)
        co_=co_+ 1;
        prev_cd = new_cd;
        iterit = randperm(N,ceil(frac_to_iter*N));
        for j=1:ceil(frac_to_iter*N)
            fun = @(r)beta*abs(I_q_theta_thetaSNew( S, K, r, theta_old , X(iterit(j),:), clst_rej, clst_acc)) + alpha * sum(KLDiv(P_h_given_x(iterit(j),:),r'));
            options = optimset('FunValCheck','on');
            %disp('Before Opt q');
            [q(:,iterit(j)),fval_new,~,~]=fminsearch(fun,q(:,iterit(j)),options);
            q(:,iterit(j)) = q(:,iterit(j))/sum(q(:,iterit(j)));
            %disp(q(:,j));
            fval_avg=fval_avg+fval_new;
            fval_latest = fval_new;
        end
        new_cd = fval_latest;
    end
    epsilon_val_final = epsilon_val 
    q = q';
   elseif nargin == 4
       q = P_h_given_x; %NXK
   else
       disp('Wrong number of arguments in EM');
   end
   %the M Step
   % Update steps as a normal GMM
   N_ks = sum(q,1);%1XK
   priors = N_ks/N; %1XK
   mu = q'*X; %KXD
   for k = 1:K mu(k,:) = mu(k,:)/N_ks(1,k); end %KXD
   for i = 1:K
       co_var_mat{i} = zeros(D,D);
       for j = 1:N
         co_var_mat{i} = co_var_mat{i} + q(j,i)*(X(j,:)-mu(i,:)).'*(X(j,:) - mu(i,:));
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
   if iter > 2
      disp(BIC(iter-1) - BIC(iter));
      if abs(BIC(iter-1) - BIC(iter)) < .001
          disp('EM has converged');
          break;
      end
   end
   iter = iter + 1;
   disp(BIC(iter-1));
end
score = BIC(iter - 1);
S = S + 1;
theta_new = {S,K,3};
% returning thr parameters of the model in theta which is a cell of size
% {K,3}
%append the now generated theta at the end of theta_old
if nargin == 7
    for s = 1:S-1
       for k = 1:K 
         theta_new{s,1,k} = theta_old{s,1,k};
         theta_new{s,2,k} = theta_old{s,2,k};
         theta_new{s,3,k} = theta_old{s,3,k};
       end  
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
