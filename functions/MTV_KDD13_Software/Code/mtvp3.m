% Function [X,iter] = mtvp3(Y,lam,rho,maxiter,tol,display)
%     total variation minimization for 3D case
%     the parallel version
%
%% Problem
%
% min 1/2||X-Y||^2 + lam\sum\sum(|x_{i,j,k}-x_{i,j+1,k}| + |x_{i,j,k}-x_{i+1,j,k}|
%     + |x_{i,j,k}-x_{i,j,k+1}|)
%
%% Input paramters:
% 
% Y-                 Observed 3D images (M x N x K)
% lam-               reguarization parameter
% rho-               the dual update length for ADMM
% maxiter-           maximum number of iteration of ADMM
% tol-               the absolute tolerance and relative tolerance
% display-           display the primal & dual error
%
%% Output paramters:
% X -                the estimated 3D image
% iter-              number of iterations
%
% For any problem, please contact with Sen Yang via senyang@asu.edu
% 
% Last modified on Aug 30, 2013
%
%% Related paper
%
% Yang, S., Wang, J., Fan, W., Zhang, X., Wonka, P., & Ye, J. An efficient
% ADMM algorithm for multidimensional anisotropic total variation 
% regularization problems. KDD (pp. 641-649). ACM, 2013