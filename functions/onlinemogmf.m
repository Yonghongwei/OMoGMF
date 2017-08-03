function [nmodel,v,label,R] =onlinemogmf(model,x,iter,w_x)
% OMoGMF main function
if (~isfield(model,'ro'))
   model.ro=1;
end
if nargin<3
   iter=3;
end
if nargin<4
   w_x=ones(size(x)) ;
end

[nmodel,label,R,h,v]=update_MoGparameters(x,model,w_x,iter);% updata MoG parmeters
nmodel=update_subspace(x,h,v,nmodel);% updata subspace
