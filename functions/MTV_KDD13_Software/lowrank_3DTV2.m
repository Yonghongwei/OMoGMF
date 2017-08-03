function [L,S]=lowrank_3DTV2(Y,sigma,lamda)
[m,n,c]=size(Y);S=zeros(m,n,c);Z=zeros(m,n,c);
DISPLAY_SETTING = 0; % display the primal & dual error 
MAX_ITER = 100;     % the maximum iteration number
TOL_SETTING = [1e-3, 1e-3]; % the absolute tolerance and relative tolerance
rho = 10;            % the dual update length for ADMM   
iter=30;mu=0.01;
for it=1:iter
[U,d,V]=svd(reshape(Y-S,[m*n,c]),'econ');
dd=diag(d);svp=sum(dd>(sigma/mu));
dd=max(dd-sigma/mu,0);
LL=U(:,1:svp)*diag(dd(1:svp))*V(:,1:svp)';
L=reshape(LL,[m,n,c]);
[S]=mtvp3(Y-L,sigma*lamda/mu,rho,MAX_ITER,TOL_SETTING,DISPLAY_SETTING);
Z=Z+mu*(Y-S-L);
mu=mu*1.3;
end

end

function Y=softhre(X,lamda)
Y=max(abs(X)-lamda,0).*sign(X);
end
