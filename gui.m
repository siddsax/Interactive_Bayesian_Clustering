function [cnt,rejected] = gui(matrix,X,K)

values=zeros(5,K);

for i=1:K
[~,sortIndex]=sort(matrix(:,K),'descend');
values(:,K)=sortIndex(1:5);
end
% will  have to check the map thing 
for i = 1:K
for j=1:5
  subplot(i,j,1),subimage(X(value(i,j),:),map);
end
end 

rejected=input('enter the rejected cluster numbers in vector form. eg [1 2 3 4 5] :','s');

