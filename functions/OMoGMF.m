function[L,E,F,label,model]= OMoGMF(model,X,iter,w_X)
%%online mog matrix factorization for video background subtraction
%input: X the data matrix;
%       model all variables and parameters for OMoGMF model 
%       model.mu the MoG pameters mu
%       model.Sigma the MoG pameters sigma^2
%       model.weight the MoG pameters pi
%       model.N controls the speed of updating MoG parameters
%       model.ro controls the speed of updating U
%       model.A and model.B auxiliary variable of subspace U
%       model.U  subspace
%       model.imgsize the frame size of video
%       model.tv.mod  switch of TV threshold: if tv.mod==1 using TV threshold
%       model.tv.lamda  parameters of TV threshold
%       iter the number of iterations
%       w_X  indicator matrix of data, if W_X_{ij}=1, X_{ij} is observed
%output:L background matrix 
%       E  residual matrix 
%       F  foreground matrix  by tv threshold
%       label data label of Gaussians 
%       model all variables and parameters for OMoGMF model 

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

if (~isfield(model,'tv'))
       model.tv.mod=0;
       model.tv.lamda=1;
end
if (~isfield(model.tv,'lamda'))
  model.tv.lamda=1;
end
if nargin<4
  w_X=ones(size(X));
end
if nargin<3
  iter=3;
end
imgsize=model.imgsize;
% lamda=tv.lamda;
[m,n]=size(X);
E=zeros(m,n);label=zeros(m,n);F=zeros(m,n);
tic
%% main outerloop
 for i=1:size(X,2)   
   if mod(i,100)==0||i==1
      disp(['Calculating the model of the ',num2str(i),'th frame']);
   end
[model,v,label0,~] =onlinemogmf(model,X(:,i),iter,w_X(:,i));% main function
 label(:,i)=label0;
 U=[model.U];
 E(:,i)=X(:,i)-U*v;
 F0=E(:,i).*(label0==1);
 %% TV threshood
 if model.tv.mod==1
FF=reshape(E(:,i),imgsize);
lam = model.tv.lamda*(model.Sigma(1));
F0= TVthre(FF,lam);
 end
 %%
F(:,i)=F0(:);
 end
%% output
L=X-E;

end
 
function F0= TVthre(FF,lam)
 if size(FF,3)==1
[F0] = mtv2(FF,lam,10,50,[1e-3, 1e-3],0);
 end
if  size(FF,3)==3
   [F0(:,:,1)]=mtvp2(FF(:,:,1), lam,10,50,[1e-3, 1e-3],0); 
   [F0(:,:,2)]=mtvp2(FF(:,:,2), lam,10,50,[1e-3, 1e-3],0); 
   [F0(:,:,3)]=mtvp2(FF(:,:,3), lam,10,50,[1e-3, 1e-3],0); 
end
end
