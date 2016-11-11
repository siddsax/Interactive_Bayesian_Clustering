% we find the probability of a single data point being generated from a
% cluster. I.E P(h|x_i) is proportional to P(x_i|h)*P(h). Here P(x_i|h) is calculaated using the normal 
% probab dist function whereas the P(h) is the prior_k found using parameters in last iteration. 
function E = P_h_givn_x(X,k,prior,mu,co_var_mat)
[n,d] = size(X);
a = (2*pi)^(0.5*d);
S = zeros(1,k);
iV = zeros(d,d,k);
for j=1:k
    if co_var_mat{j} == zeros(d,d), co_var_mat{j}=ones(d,d)*eps;  end
    S(j) = sqrt(det(co_var_mat{j}));
    iV(:,:,j) = inv(co_var_mat{j});    
end
E = zeros(n,k);
for i=1:n    
    for j=1:k
        dXM = X(i,:)-mu(j,:);
        pl = exp(-0.5*dXM*iV(:,:,j)*dXM')/(a*S(j));
        E(i,j) = prior(1,j)*pl;
    end
    E(i,:) = E(i,:)/sum(E(i,:));
end


