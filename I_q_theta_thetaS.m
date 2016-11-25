% this finction gives back the I_q for a particular s as mentioned in the paper with the input
% parameters as the the the reject-accept iteration s, the number of
% clusters, the arbitrary probab, the co variance matrices of all the last
% steps for all clusters, the priors (prior is pi_k as mentioned in our class, it is hence a 
% 1XK array) of all the last steps for all clusters, data X, the indices of rejected and 
% accepted clusters)
function [ req_p ] = I_q_theta_thetaS( s, K, q, co_var_mat_s, prior_s, mu_s, X, clust_rej, clust_acc )
req_p = 0;
[N,~] = size(X);
P_h_given_x = P_h_givn_x(X,K,prior_s,mu_s,co_var_mat_s); %NXK
if isreal(P_h_given_x) == 0
    disp('here line 12');
    pause;
end    
P_q_hs = P_h_hs(q,P_h_given_x,N); %KXK
if isreal(P_q_hs) == 0
    disp('here line 17');
    pause;
end
Marg_P_h = sum(P_h_given_x,1)/N; %1XK
for i = 1:K
    marg_prob_q = 0;
    for j = 1:N
        marg_prob_q = marg_prob_q + q(i,j); 
    end    
    marg_prob_q = marg_prob_q/N;
    if ~isempty(clust_acc)
       for j = 1:size(clust_acc(1,:),2)
          req_p = req_p - P_q_hs(i,clust_acc(1,j))*(log(P_q_hs(i,clust_acc(1,j))) - log(Marg_P_h(1,clust_acc(1,j))) - log(marg_prob_q));
          if isreal(req_p) == 0
              disp('here line 31');
              pause;
          end
       end   
    end   
    if ~isempty(clust_rej) 
        for j = 1:size(clust_rej(1,:),2)
            req_p = req_p + P_q_hs(i,clust_rej(1,j))*(log(P_q_hs(i,clust_rej(1,j))) - log(Marg_P_h(1,clust_rej(1,j)))-log(marg_prob_q));
            if isreal(req_p) == 0
               disp('here line 39');
               pause;
            end
        end
    end   
    if isreal(req_p) == 0
        disp('here line 45');
        pause;
    end    
end

