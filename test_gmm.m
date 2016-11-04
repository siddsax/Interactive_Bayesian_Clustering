close all;
X = load('test_kmeans.txt');

fprintf('First, we just plot some data...\n\n');

figure(1); title('figure 1');
plot_gmm(X);

mypause;

fprintf(['Now, we test the K=2 implementation of gmm; if this looks\n' ...
         'totally unreasonable, you probably have a bug.  We use (totally)\n' ...
         'random initialization...\n\n']);


figure(2); title('figure 2');
[mu,pk,z,si2,CLL,ILL,BIC] = gmm(X,2);
plot_gmm(X,mu,z,BIC);

mypause;

fprintf(['Now, we run sixteen times with different initializations and\n' ...
         'compare the BIC (lower BIC is better).  For this, we will run with K=3, which seems\n' ...
         'more reasonable for this data...\n\n']);

figure(3); title('figure 3');
for ii=1:16,
  [mu,pk,z,si2,CLL,ILL,BIC] = gmm(X,3);
  subplot(4,4,ii);
  plot_gmm(X,mu,z,BIC);
end;
mypause;

fprintf('Finally, we consider BIC as a function of K (best K should have smallest BIC)...\n\n');

allK = [2 3 4 5 6 8 10 15 20];
allS = Inf + zeros(size(allK));
allMu = {}; allZ = {};
for ii=1:length(allK),
  fprintf('  k=%d...', allK(ii));
  for rep=1:20,
    [mu,pk,z,si2,CLL,ILL,BIC] = gmm(X,allK(ii));
    if BIC < allS(ii),
      allS(ii) = BIC;
      allMu{ii} = mu;
      allZ{ii} = z;
    end;
  end;
end;
fprintf('\n');
figure(4); title('figure 4');
plot(allK,allS,'bx-');
xlabel('K'); ylabel('BIC');


figure(5); title('figure 5');
for ii=1:5,
  subplot(2,3,ii);
  plot_gmm(X,allMu{ii},allZ{ii},allS(ii));
end;
mypause;
