function [wo,epsreturn] = calculateweightgeneral(T,d,lvec,dp)

%alternate form - tries to minimize epsilon

%%%%%%compute the optimal weight%%%%%%%%%%%%%%
%kvec=sk:k;
N=floor(T/2);
M=T-N;
kvec = floor(lvec*sqrt(T)/2);

%optimal weight - kks
%coeff = [0 1 2 3 4 5 6 7 8 9 10];
%x = [1 0 0 0 0 0 0 0 0 0 0]';

coeff = (2:(floor(dp)+2))-1;
%coeff = (1:(d/2+1))-1;
x = zeros(size(coeff));x(1) = 1;

kmat(1,:) = ones(size(kvec));
sc = 1*(1/T)^(1/2);
Aeq(1,:) = kmat(1,:);
Aeq(1,end+1) = 0;
Beq = 1;
count = 0;
options=optimset('Display','off');
for i=2:length(x)
    %kmat(i,:) = ((gamma(kvec+1-al+coeff(i)/d)./gamma(kvec))) * (gamma(T)/gamma(T+coeff(2)/d));
    kmat(i,:) = (kvec).^(coeff(i)/d) * (M^(-coeff(i)/d));
    count=count+1;
    A(count,1:length(kvec)) = kmat(i,:);
    A(count,length(kvec)+1) = -1;
    B(count) = 0;
    count=count+1;
    A(count,1:length(kvec)) = -kmat(i,:);
    A(count,length(kvec)+1) = -1;
    B(count) = 0;
end
for i=0:d
   A0(:,i+1) = lvec.^i; 
end
A1 = A0(:,2:end);
eta1 = sqrt(det(A1'*A1)/det(A0'*A0));

eta2=(1/((1.5)^(1/d)-1));

eta = eta2%max(eta1,eta2) % ensures large enough eta for small d

%10*sqrt(length(kvec));%*T^(0.4);
winit = ones(size(kvec))/length(kvec);
wleft = winit;wleft(end+1) = -100;0;
wright = winit;wright(end+1) = 100;%((k/M)^(2/d))/kap;
winit(end+1) = sc;
wo = fmincon(@(ws) ws(end), winit, A, B', Aeq, Beq,-eta*abs(wleft),eta*abs(wright),[],options);
epsreturn = wo(end);
wo(end)=[];

