%simple compare

clc
clear all

matlabpool close force local

matlabpool(7)

matlabpool size

p=0.65; d=6; a=2.25;b=2.15;
q=1-p;

%-------f density function------------%
f = @(x) p*prod(betapdf((x),a,b),2)+q;
%-------------------------------------%

T=10^6;
P=sum(rand(T,1)<p);
Q=T-P;
Xb=betarnd(a,b,P,d);
Xu=rand(Q,d);
X=[Xb;Xu];

X=X(randperm(T),:);

%e=mean(-log(f(X)));
alpha=0.5;
e=mean((f(X)).^(alpha-1));

nn=8;

Tvec = floor(logspace(log10(1000),log10(10000),nn));
epsvector=zeros(nn,1);



for tloop = 1:nn

    tloop

Niter=28;
%tic
lenofvec=20;lenofvecmid=round(lenofvec/2);
T=Tvec(tloop);dp=d;fraction = 0.5;
kmax = 2;%max(5,100*floor((T/2)^(1/(1+d))))
avec = linspace(.3,3,lenofvec);
%kvec = floor(avec*sqrt(T/2));
kvec = floor(avec*((T/2)^(2/(2+2))));
[wo,epsval] = calculateweightgeneral(T,d,avec,dp);
epsvector(tloop) = epsval;errorvec=[0,0,0,0,0,0];errorvecb=[0,0,0,0,0,0];
kvl = length(kvec);
errorvecT = zeros(Niter,6);
parfor iter=1:Niter
    
    %iter
    
P=sum(rand(T,1)<p);Q=T-P;Xb=betarnd(a,b,P,d);
Xu=rand(Q,d);X=[Xb;Xu];X=X(randperm(T),:);

G_plug = zeros(kvl,1);
G_BC = zeros(kvl,1);

for i=1:kvl
[g_plug,g_BC] = truncatedRenyiestimate(X,kvec(i),alpha);
%[g_plug,g_BC] = truncatedShannnonestimate(X,kvec(i));
G_plug(i) = g_plug;
G_BC(i) = g_BC;
end

G_w = wo*G_plug;
G_wa = wo*G_BC;

[G_plug,G_BC] = truncatedRenyiestimate(X,kmax,alpha);
%G_plug = G_plug(end);
%G_BC = G_BC(end);

[G_plug_knn,G_BC_knn,G_h_knn] = knn_Renyi_entropy_estimate(X,2,kmax,fraction,alpha);

errorvecT(iter,:) = (([G_plug,G_BC,G_plug_knn,G_BC_knn,G_h_knn,G_wa]-e));
errorvec = errorvec+(([G_plug,G_BC,G_plug_knn,G_BC_knn,G_h_knn,G_wa]-e).^2);
errorvecb = errorvecb+(([G_plug,G_BC,G_plug_knn,G_BC_knn,G_h_knn,G_wa]-e));
end
errorvec = errorvec/Niter;
errorvecb=errorvecb/Niter;

%toc

    errorTvec(tloop,:)=errorvec
    errorTvecStd(tloop,:)=sqrt(var(errorvecT))
    errorTvecb(tloop,:)=errorvecb

end

matlabpool close 

log10(errorTvec(1,1))-log10(errorTvec(end,1))
log10(errorTvec(1,2))-log10(errorTvec(end,2))
log10(errorTvec(1,3))-log10(errorTvec(end,3))
log10(errorTvec(1,4))-log10(errorTvec(end,4))
log10(errorTvec(1,5))-log10(errorTvec(end,5))
log10(errorTvec(1,6))-log10(errorTvec(end,6))



