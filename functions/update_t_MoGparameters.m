function [model,label,R,h,v,tx,tau]=update_t_MoGparameters(x,model,tau,iter)
U=model.U;
N=model.N;
nmodel.weight=model.weight;
nmodel.mu=model.mu;
nmodel.Sigma=model.Sigma;
r=size(U,2);

m=size(x,1);
h=ones(m,1);
if nargin<4
iter=20;
end
imgsize=model.imgsize;
if (~isfield(model,'v'))
v=(U'.*repmat(h',r,1)*U+0.00001*eye(r))^-1*U'*(h.*x);
else
v=model.v;
end

if model.preAlignment==1;
[B,tau] = regMGNC(reshape(U*v,imgsize),reshape(x,imgsize),tau,fix(log(m)/log(2)/2-3),2);
else
[B,tau] = regImg(reshape(U*v,imgsize),reshape(x,imgsize),tau,reshape(h,imgsize),1);
end

tx=B(:);d=size(tau,1);
%% main loop
for tt=1:iter
J =getJ(reshape(x,imgsize),tau);
 dtau=zeros(size(tau));T=[U,-J];
  for itt=1:1
%% E_step
  [R] = expectation((tx-U*v)', nmodel);
  [~,label]=max(R,[],2);
  label=label';
%% M_step
%M_step for w sigma N
  [nmodel] = maximizationModel((tx-U*v)',R,model,N);
%M_step for v dtau
  h=sum(R.*repmat(1./(2*nmodel.Sigma),m,1),2);
  temp=T'.*repmat(h',[r+d,1]);
  temp2=(temp*T+0.00001*eye(r+d))\temp*tx;
  v=temp2(1:r);dtau=temp2((r+1):end);
  end
%% updata tau
tau=tau+dtau;
[B] = warpImg(reshape(x,imgsize),tau);tx=B(:);
if sum(abs(dtau))<0.0001
break;
end
end
model.weight=nmodel.weight;
model.mu=nmodel.mu;
model.Sigma=nmodel.Sigma;
model.v=v;

end

function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;
end

function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
logRho = zeros(n,k);


for i = 1:k
    logRho(:,i) = loggausspdf(X,mu(i),Sigma(i));
end
logRho = bsxfun(@plus,logRho,log(w));
T = logsumexp(logRho,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);
end

function [nmodel] = maximizationModel(X,R,model,N)
ww=model.weight;
mmu=model.mu;
ssigma=model.Sigma;
% [d,n] = size(X);
k = size(R,2);
nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
% mu = bsxfun(@times, X*R, 1./nk);
w = nk/size(R,1);
Sigma = zeros(1,k);
sqrtR = sqrt(R);
new_Nk=ww*N+sum(R);
new_N=N+size(X,2);
new_w=new_Nk/new_N;
for i = 1:k
    Xo = bsxfun(@minus,X,mu(i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(i) = (Xo*Xo'+ssigma(i)*ww(i)*N)/new_Nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end
nmodel.mu = mu;
nmodel.Sigma =Sigma;
nmodel.weight = new_w;
end



