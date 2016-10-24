% we find the probability of a single data point being generated from a
% cluster. I.E P(h|x_i) is proportional to P(x_i|h)*P(h). Here P(x_i|h) is calculaated using the normal 
% probab dist function whereas the P(h) is the prior_k found using parameters in last iteration. 
function [ req_prob ] = P_h_givn_x( dat_pt_indice, X, co_var_mat_k, prior_k, mu_k )
[~,D] = size(X); 
part1 = power(1/(2*pi),D/2.0)*power(Det(co_var_mat_k),1/2.0);
part2 = exp((-1/2)*((X(dat_pt_indice,:) - mu_k).')*(co_var_mat_k\(X(dat_pt_indice,:) - mu_k)));
normal_part = (part1.^-1)*part2;
req_prob = normal_part*prior_k;
end


