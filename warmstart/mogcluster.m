function [model label,R] =mogcluster(x,k)
R = initialization(x',k);%这里有错吧
[~,label(1,:)] = max(R,[],2);%得到第n个是哪一类的。
R = R(:,unique(label));
model.mu = zeros(1,k);
model.Sigma =mean((x).^2)*k*10*0.1.^(1:k);
model.Sigma=-sort(-model.Sigma);
model.weight = ones(1,k)/k;
% nk = sum(R,1);
% model.weight = nk/size(R,1);
for tt=1:50
[R, llh(tt)] = expectation(x', model);
[model] = maximizationModel(x',R);
[~,label]=max(R');
% label=label';
%stop 
tol=0.000001;
 if tt>1
if abs(llh(tt)-llh(tt-1)) < tol*abs(llh(tt))
    fprintf('Converged in %d steps.\n',tt-1);
     break;
 end
end
end

end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    [u,~,label] = unique(label);
    while k ~= length(u)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
        [u,~,label] = unique(label);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d  %initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end
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


function [model] = maximizationModel(X,R)

[d,n] = size(X);
k = size(R,2);

nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
% mu = bsxfun(@times, X*R, 1./nk);
w = nk/size(R,1);
Sigma = zeros(1,k);
sqrtR = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(i) = Xo*Xo'/nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;

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

