function [G_plug,G_BC,G_u,G_w] = weighted_Renyi_entropy_estimate(S,sk,k,kvec,wo,fraction,alpha)

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
%G_u: uniformly weighted and evaluated over range kvec
%G_w: optimally weight and evaluated over range kvec


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



%----------------------------------
k2=max(kvec);N=floor(T/2);
M=T-N;
%weighted estimators
    Sn = S(1:N,:);
    Sm = S((N+1):T,:);
    %obtain kNN distances
    [nnidx,dists]=annquery(Sm',Sn',k2);
    kdist=dists(kvec,:)';
    
    %note gamma(kvec)./gamma(kvec+1-alpha) = beta(kvec,1-alpha)./gamma(1-alpha)
    etahatvec = (beta(kvec,1-alpha)./gamma(1-alpha)).*mean((((M).*cdunitball.*(kdist.^d))).^(1-alpha));
    
    %uniform estimator
    wu = ones(size(wo))/length(wo);
    G_u = (sum(wu.*etahatvec));
    
    %optimally weighted estimator
    G_w = (sum(wo.*etahatvec));
  
    
