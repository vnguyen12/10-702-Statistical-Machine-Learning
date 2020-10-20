function [G_plug,G_BC,G_h] = knn_Renyi_entropy_estimate(S,sk,k,fraction,alpha)

%Function computes \int g(f) f d\mu entropy for different choice of
%estimators with g(u) = u^{\alpha-1}.

%Input: 
%(i) S is T(sample size) x d(intrinsic dimension) data matrix
%sk: number of neighbors k for estimating G_BC
%k: number of neighbors for estimating G_plug
%kvec: vector of k values for estimating G_u and G_w
%wo: optimal weight computed using calculateweightgeneral or calculateweightBC.m
%fraction: fraction of N/(N+M) for the plug-in estimator
%alpha: value of alpha

%Output: 
%G_plug: Single weight estimator (Plug-in estimator) (evaluated at k)
%G_BC: Single weight with bias correction (Leonenko 2008 estimator) (evaluated at sk)
%G_h: entropic graph estimator of Hero, Pal (evaluated at sk)


[T,d] = size(S);

%Volume of unit ball in d dimensions
cdunitball=(pi^(d/2))/gamma(d/2+1);

N=floor(T*fraction);
M=T-N;

Sn = S(1:N,:);
Sm = S((N+1):T,:);

%obtain kNN distances
K=max(sk,k);
[nnidx,dists]=annquery(Sm',Sn',K);


%Plug-in estimator
kdist=dists(k,:)';
fhat_inv = (1./(k-1)).*((M).*cdunitball.*kdist.^d);
G_plug = (mean((fhat_inv).^(1-alpha)));

%Bias corrected estimator
kdist=dists(sk,:)';
fhat_inv = (1./(sk-1)).*((M).*cdunitball.*kdist.^d);
G_BC = (gamma(sk)/gamma(sk+1-alpha)).*((sk-1).^(1-alpha)).*mean((fhat_inv).^(1-alpha));


%
gam = d*(1-alpha);
[nnidx,dists]=annquery(S',S',sk+1);
dists = dists(1:sk,:);
Lf = sum(sum(dists.^(gam)));

U = rand(T,d);
[nnidx,dists]=annquery(U',U',sk+1);
dists = dists(1:sk,:);
Lunif = sum(sum(dists.^(gam)));

G_h = Lf/Lunif;
