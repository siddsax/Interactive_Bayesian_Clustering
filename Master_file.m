%file_name_train = input('Please mention the path to pick up clustering data');
%K = input('The desired number of clusters');
%eps = input('Tolerance of the data for EM');
eps = 10^-20;
%max_iter = input('maximum iterations for EM');
max_iter = 1000;
%load('C:\Users\sidds\OneDrive\Documents\MATLAB\Project\cifar-10-batches-mat\data_batch_1.mat');
rng default; % For reproducibility
X = [randn(100,2)*0.75+ones(100,2);
    randn(100,2)*0.5-ones(100,2)];
figure;
plot(X(:,1),X(:,2),'.');
title 'Randomly Generated Data';
[N,D] = size(X);
K = 3;
%%
[theta,BIC] = EM(X, max_iter, K,eps);
S = 1;
co_var_mat = {K};
mu = zeros(K,D);
prior = zeros(1,K);
for k = 1:K
   co_var_mat{k} = theta{S,1,k};
   mu(k,:) = theta{S,2,k};
   prior(1,k) = theta{S,3,k};
end
Probabilities = P_h_givn_x(X, K,prior, mu, co_var_mat);
figure(2); title('figure 2');
plot_gmm(X,mu,Probabilities,BIC);
S = S + 1;
%%
 acc_clst = [];
 rej_clst = [];
 str = input('Do you want feedback Y/N','s');
 if eq(str,'Y')
     for z = 1:K
         indice = input('accepted cluster number, press 0 when you are done');
         if indice == 0 break; end
         acc_clst = [acc_clst,indice];
     end
     for z = 1:K
         indice = input('rejected cluster number, press 0 when you are done');
         if indice == 0 break; end
         rej_clst = [rej_clst,indice];
     end
 end    
 [theta,BIC] = EM(X, max_iter, K,eps,theta, rej_clst, acc_clst);
 for k = 1:K
   co_var_mat{k} = theta{S,1,k};
   mu(k,:) = theta{S,2,k};
   prior(1,k) = theta{S,3,k};
end
Probabilities = P_h_givn_x(X, K,prior, mu, co_var_mat);

figure(5); title('figure 5');
plot_gmm(X,mu,Probabilities,BIC);
S = S + 1;