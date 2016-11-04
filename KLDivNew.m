function [ sum_KL ] = KLDivNew(P,Q,K,N)
sum_KL=0;
for i = 1:N
  for j = 1:K
      sum_KL = sum_KL + KLDiv(P(j,i),Q(j,i));
  end
end
end
