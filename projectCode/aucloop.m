
% code to generate Area under ROC curves
clc
clear all

matlabpool close force local

matlabpool(7)

matlabpool size

Titer = 500; %total number of points

zerovec = zeros(Titer/2,1);
onesvec = ones(Titer/2,1);

zeon = [zerovec;onesvec];
zeon = zeon(randperm(Titer)');


deltavec = .1:.1:1;

for deltai = 1:length(deltavec)


T=floor(50*10^2);d=6;a=10;b=10;dp=d/2;a2=10-deltavec(deltai);b2=10-deltavec(deltai);p=0.75;
sk = floor(0.3*sqrt(T));k=floor(2*sqrt(T));
[wo,epsval] = calculateweightshannon_lowervar(T,d,sk,k,dp);


shannonsimp = zeros(Titer,1);shannonsi=shannonsimp;
shannonweight = shannonsimp;
maxdist=zeros(Titer,1);
mstcount=zeros(Titer,1);
spacingcount=zeros(Titer,1);

%true entropy
q=1-p;
%-------f density function------------%
f = @(x) p*prod(betapdf((x),a,b),2)+q;
%-------------------------------------%
%-------f density function------------%
f2 = @(x) p*prod(betapdf((x),a2,b2),2)+q;
%-------------------------------------%


e = shannonsimp;ll=shannonsimp;
    TT=10^6;
    P=sum(rand(TT,1)<p);Q=TT-P;
        Xb=betarnd(a,b,P,d);Xu=rand(Q,d);X=[Xb;Xu];X=X(randperm(TT),:);
        e1 = -mean(log(f(X)));
        P=sum(rand(TT,1)<p);Q=TT-P;
        Xb=betarnd(a2,b2,P,d);Xu=rand(Q,d);X=[Xb;Xu];X=X(randperm(TT),:);
        e2 = -mean(log(f2(X)));

parfor i = 1:Titer
   %i
    if(zeon(i)==0)
        P=sum(rand(T,1)<p);Q=T-P;
        Xb=betarnd(a,b,P,d);Xu=rand(Q,d);X=[Xb;Xu];X=X(randperm(T),:);
        [shannonsi(i),shannonsimp(i),shannonweight(i),maxdist(i),mstcount(i),spacingcount(i)] = testforuniformity(X,wo);
        e(i) = e1;
        ll(i) = sum(log(f2(X))) - sum(log(f(X)));
    else
        P=sum(rand(T,1)<p);Q=T-P;
        Xb=betarnd(a2,b2,P,d);Xu=rand(Q,d);X=[Xb;Xu];X=X(randperm(T),:);
        [shannonsi(i),shannonsimp(i),shannonweight(i),maxdist(i),mstcount(i),spacingcount(i)] = testforuniformity(X,wo);
        e(i)=e2;
        ll(i) = sum(log(f2(X))) - sum(log(f(X)));
    end
    
end
%ll=ll/max(ll);


%%

[Xk,Yk] = perfcurve(zeon,ll,1); np(deltai) = areaundercurve(Xk',Yk')

[Xk,Yk] = perfcurve(zeon,shannonsi,1); ssi(deltai) = areaundercurve(Xk',Yk')

[Xk,Yk] = perfcurve(zeon,shannonsimp,1); ssism(deltai) = areaundercurve(Xk',Yk')


[Xk,Yk] = perfcurve(zeon,shannonweight,1); sswgt(deltai) = areaundercurve(Xk',Yk')

end
matlabpool close 
%%
npS = smooth(np);
ssiS = smooth(ssi);%ssiS(end) = 0.92;
ssismS=smooth(ssism);%ssismS(end)=0.94;
sswgtS=smooth(sswgt);%sswgtS(end) = 0.96;
plot(deltavec,npS,'k',deltavec,ssiS,'r',deltavec,ssismS,'b',deltavec,sswgtS,'g')
legend('Neyman-Pearson test','Standard kernel plug-in estimate','Truncated kernel plug-in estimate','Weighted estimate')
xlabel('$\delta$')
ylabel('Area under ROC curve')



















