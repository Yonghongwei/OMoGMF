function [U,V]=initi_OMoG(XX,r)
X=reshape(permute(XX,[1,3,2]),[size(XX,1)*size(XX,3),size(XX,2)]);
[U,V]=PCA(X,r);beta=0.0001;Y=zeros(size(X));E=X-U*V';
ro=1.5;iter=50;mu=0.01;
for i=1:iter
V=((X-E+Y/mu)'*U)/(U'*U+beta/mu*eye(r));
U=(X-E+Y/mu)*V/(V'*V+beta/mu*eye(r));
E=threL1(X-U*V'+Y/mu,1/mu);
Y=Y+mu*(X-U*V'-E);
mu=mu*ro;
end
end
function Y=threL1(X,lamda)
Y=(abs(X)-lamda).*sign(X);
end