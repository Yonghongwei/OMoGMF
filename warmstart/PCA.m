function [Y,W,d,flag] = PCA(X,k)
% PCA computes the k principal components with largest associated
% eigenvalue.
%
% INPUTS
% X: (n x d) zero mean data matrix (samples are rows)
% k: number of principal components
%
% OUTPUTS
% Y: (n x k) first k principal components
% W: (d x k) respective loadings of PCs
% d: eigenvalues in descending order
% flag: convergence flag (see eigs)
%
% (c) 2008 Christian Sigg (chrsigg@inf.ethz.ch)
% See license file for details.

% eigs options
options.issym = 1;
options.disp = 0;

r = rank(X);
if nargin < 2
    k = r;
end

if(size(X,1)>=size(X,2))
    if (k == r)
        [W,D] = eig(X'*X);
        flag = 1;
    else
        [W,D,flag] = eigs(X'*X,k,'lm',options);
    end
    [d,indx] = sort(diag(D),'descend');
    Y = X*W;
else
    if (k == r)
        [W,D] = eig(X*X');
        flag = 1;        
    else
        [W,D,flag] = eigs(X*X',k,'lm',options);
    end
    [d,indx] = sort(diag(D),'descend');
    W = X'*W;
    for j=1:size(W,2)
	W(:,j) = W(:,j)/norm(W(:,j));
    end
    Y = X*W;
end

Y = Y(:,indx(1:k));
W = W(:,indx(1:k));
d = d(1:k);
end
