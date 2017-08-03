function [model] = warmstart(X,r,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');
 disp('initializing U...');
%% PCA and L1MF for warmstart
[U,V]=initi_OMoG(X,r);
%% Calculating the model
disp('initializing MoG parameters...');
L=U*V';noi=X-L;
% [nmodel,~,R] =mogcluster(noi,k); %using GMM for  initialization of mog parameters
nmodel.mu = zeros(1,k);
nmodel.Sigma =var(noi(:))*50*0.05.^(1:k); %a simple initialization of mog parameters
nmodel.Sigma=-sort(-nmodel.Sigma);
nmodel.weight = ones(1,k)/k;
%% Calculating A and B
W=1./(abs(noi)+0.001);
m=size(X,1);
A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(V'*diag(W(i,:))*V+0.001*eye(r))^-1*0.02;
    B(:,i)=V'*diag(W(i,:))*(X(i,:))'/0.02;
end
disp('Warm start is over. ');
%% output
model.Sigma=nmodel.Sigma;
model.weight=nmodel.weight;
model.mu=nmodel.mu;
model.A=A;
model.B=B;
model.U=U;
end