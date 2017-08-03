function [model] = t_warmstart(X,r,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');

%% prealignment
 disp('prealignment...');
[m,n,c]=size(X);
TX=X;
for i=1:3
    if i==1
[TX,Tau] = preAlign(X,median(TX,3));
    else
[TX,Tau] = preAlign(X,median(TX,3),Tau);
    end
end
%% initialization
TX=reshape(TX,[m*n,c]);
X=reshape(TX,[m*n,c]);
[U,V]=initi_OMoG(TX,r);
noise=TX-U*V';
W=1./(abs(noise)+0.001);
A=zeros(r,r,m*n);B=zeros(r,m*n);
for i=1:m*n
    A(:,:,i)=(V'*diag(W(i,:))*V+0.001*eye(r))^-1*0.1;
    B(:,i)=V'*diag(W(i,:))*(TX(i,:))'/0.1;
end

model.Sigma=(var(noise(:)))*50*0.05.^(1:k);
model.Sigma=-sort(-model.Sigma);
model.weight=ones(1,k)/k;
model.mu=zeros(1,k);

model.A=A;
model.B=B;
model.U=U;
model.N=20;model.ro=0.99;
model.imgsize=[m,n];model.preAlignment=0;
%Tau=zeros(6,size(X,2));
 disp('Align image...');
tic;
%% alignment
for i=1:3
[~,~,~,~,model,Tau]=  t_OMoGMF(model,X,2,Tau);
end
disp('Warm start is over. ');
end