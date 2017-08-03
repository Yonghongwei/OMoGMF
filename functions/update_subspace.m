function model=update_subspace(x,h,v,model)
%% updata auxiliary variable A and B for U
A=model.A;
B=model.B;
ro=model.ro;
%% tensor coding for upspeeding
temp= 1/ro*sum(bsxfun(@times,A,v'),2);
hh=permute(h,[3,2,1]);
new_A=1/ro*A-bsxfun(@rdivide,bsxfun(@times,bsxfun(@times,hh,temp),permute(temp,[2,1,3])),...
    (1+bsxfun(@times,hh,sum(bsxfun(@times,v,temp),1))));
new_B=B*ro+bsxfun(@times,(h.*x)',v);
U=permute(sum(bsxfun(@times,new_A,permute(new_B,[1,3,2])),1), [3,2,1]);
%% output 
model.A=new_A;
model.B=new_B;
model.U=U;
