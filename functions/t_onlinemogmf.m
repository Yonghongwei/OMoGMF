function [nmodel,v,label,R,B,tau] =t_onlinemogmf(model,x,tau,iter)

if (~isfield(model,'ro'))
   model.ro=1;
end

if nargin<3
  tau=zeros(6,1);  
end
[nmodel,label,R,h,v,B,tau]=update_t_MoGparameters(x,model,tau,iter);
% nmodel.mod=2;
nmodel=update_subspace(B(:),h,v,nmodel);
