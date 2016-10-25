%file_name_train = input('Please mention the path to pick up clustering data');
%K = input('The desired number of clusters');
K = 10;
%eps = input('Tolerance of the data for EM');
eps = .0001;
%max_iter = input('maximum iterations for EM');
max_iter = 1000;
load('C:\Users\sidds\OneDrive\Documents\MATLAB\Project\cifar-10-batches-mat\data_batch_1.mat');
T = data;
Y = labels;
X = im2double(T);
[~,theta] = EM(X, max_iter, K);
[N,~] = size(X);
Probabilities = zeros(N,K);
S = 1;
for i = 1:K
    for j = 1:N
       Probabilities(j,i) = P_h_givn_x(j, X, theta{S,i,1}, theta{S,i,3}, theta{S,i,2});
    end
end
 for i = 1:K
    c=flipud(unique(sort(Probabilities(:,i))));
    result=c(1:5); 
    ind=find(Probabilities(:,i)>=c(10));
    disp(ind);
 end    