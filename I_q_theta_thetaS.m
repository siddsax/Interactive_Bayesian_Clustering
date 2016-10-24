% this finction gives back the I_q for a particular s as mentioned in the paper with the input
% parameters as the the the reject-accept iteration s, the number of
% clusters, the arbitrary probab, the co variance matrices of all the last
% steps for all clusters, the priors (prior is pi_k as mentioned in our class) of all the last steps for all
% clusters, data X, the indices of rejected and accepted clusters)
function [ req_p ] = I_q_theta_thetaS( s, K, q,co_var_mat, prior, mu, X, clst_rej, clst_acc )
req_p = 0;
for i = 1:K
    marg_prob_q = 0;
    for j = 1:N
        marg_prob_q = marg_prob_q + q(1,j); 
    end    
    for j = 1:size(clst_acc(s,:),1)
        req_p = req_p - p_h_hs(q,co_var_mat(s,clust_acc(s,j)),prior(s,clust_acc(s,j)), X,mu(s,clust_acc(s,j)))*log(p_h_hs(q,co_var_mat(s,clust_acc(s,j)),prior(s,clust_acc(s,j)), X, mu(s,clust_acc(s,j)))/(marg_prob_h(X,co_var_mat(s,clst_acc(s,j)),prior(s,clst_acc(s,j)),mu(s,clst_acc(s,j)))*marg_prob_q));
    end
    for j = 1:size(clst_rej(s,:),1)
        req_p = req_p + p_h_hs(q,co_var_mat(s,clust_rej(s,j)),prior(s,clust_rej(s,j)), X,mu(s,clust_rej(s,j)))*log(p_h_hs(q,co_var_mat(s,clust_rej(s,j)),prior(s,clust_rej(s,j)), X, mu(s,clust_rej(s,j)))/(marg_prob_h(X,co_var_mat(s,clst_rej(s,j)),prior(s,clst_rej(s,j)),mu(s,clst_rej(s,j)))*marg_prob_q));
    end
end    
end

