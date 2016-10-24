function [ sum_KL ] = KLDivNew(P,Q)
sum_KL=0;
[K,N]=size(q);
for i = 1:N
  for j = 1:K
      sum_KL = sum_KL + KLDiv(P(j,i),Q(j,i));
  end
end
end
