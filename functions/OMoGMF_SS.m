function [L,model]= OMoGMF_SS(model,X,iter)
%%online mog matrix factorization for video background subtraction with subsampling
%%
%input:X the data matrix;
%      model
%      model.N controls the speed of updating MoG parameters
%      model.ro controls the speed of updating U
%      model.Sigma the MoG pameters sigma^2
%      model.weight the MoG pameters pi
%      model.mu the MoG pameters mu
%      model.A and model.B auxiliary variable of subspace U
%      model.U  subspace
%      model.SSrate  subsampling rate
%      model.imgsize the frame size of video
%      iter the number of iterations
%      W_X  indicator matrix of data, if W_X_{ij}=1, X_{ij} is observed
%output:L background matrix 
%       E  residual matrix 
%       F  foreground matrix  by tv threshold
%       U  subspace
%       label data label of Gaussians 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the following papers:
% [1] Deyu Meng, Fernando De la Toree, Matrix Factorization with Unknown Noise. ICCV, 2013.
% [2] Hongwei Yong, Deyu Meng, Wangmeng Zuo, Lei Zhang, Robust Online Matrix Factorization for Dynamic Background Subtraction, IEEE Transactions on Pattern Analysis  and Machine Intelligence (TPAMI), 2017. In press.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resources for research purpose only, shall not be used for commercial purposes! All copyrights belong to the original anthors. The technology has applied for patents.  If you want to purchase the patents for commercial purposes, please contact the corresponding author: Deyu Meng, dymeng@mail.xjtu.edu.cn. Thank you!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Hongwei Yong. If having any question, feel free to contact: cshyong@comp.polyu.edu.hk or yonghw@stu.xjtu.edu.cn.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.0, release date: 2017.8.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isfield(model,'ro'))
   model.ro=0.999;
end
if nargin<3
  iter=3;
end
if (~isfield(model,'SSrate'))
   model.SSrate=1;
end
L=zeros(size(X));
tic
%% main outerloop
 for i=1:size(X,2)
    Ind=randperm(size(X,1));
     ind=Ind(1:fix(model.SSrate*size(X,1)));
if mod(i,200)==0||i==1
      disp(['Calculating the model of the ',num2str(i),'th frame']);
end
nmodel=model;
nmodel.U=model.U(ind,:);
nmodel.A=model.A(:,:,ind);
nmodel.B=model.B(:,ind);
[nmodel,v] =onlinemogmf(nmodel,X(ind,i),iter);
model.U(ind,:)=nmodel.U;
model.A(:,:,ind)=nmodel.A;
model.B(:,ind)=nmodel.B;       
model.weight=nmodel.weight;
model.mu=nmodel.mu;
model.Sigma=nmodel.Sigma;        
U=model.U;
L(:,i)=U*v;
 end


