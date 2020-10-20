%simple code to illustrate trucated kernel density estimate

clc
clear all

% simple kernel 

T=10^3;
d=6;

N=T/2;
M=T-N;

%number of nearest neighbors
k=10;
%radius
dia = (k/M)^(1/d);
rad = dia/2;

%draw sample realizations
X = rand(T,d);

%partition
Xr = X(1:N,:);
Xq = X(((N+1):T), :);

%get lval

inds = zeros(N,M);

for i=1:N
    inds(i,:) = (max(abs(Xq-repmat(Xr(i,:),M,1)),[],2)<rad);
end

lval = sum(inds~=0,2); 

%now get truncated volume for Xr
radc = ((k/M)^(1/d))*ones(size(lval));
vol = ones(size(lval));
for i=1:d
    vol = vol.*((min([abs(Xr(:,i)-zeros(N,1)) abs(Xr(:,i)-ones(N,1)) rad*ones(N,1)],[],2)) + radc/2);
end

%length(find(vol==(radc.^d)))/N


%now construct truncated kde estimate

fhat_unc = lval./(M*(radc.^d));
fhat_c = lval./(M*vol);

fhat_unc(fhat_unc==0) = 1/k;
fhat_c(fhat_c==0) = 1/k;


mean(log(fhat_unc))
mean(log(fhat_c))


