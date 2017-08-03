% This is the testing demo for video background subtraction with OMoGMF
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

clear all
currentFolder = pwd;
addpath(genpath(currentFolder))
%% read data:
% X :data matrix, each coloum vector is one frame data
% imgsize: size of each frame image
load 'test_data\airport.mat'
imgsize=[size(X,1),size(X,2)];
X=im2double(reshape(X,[size(X,1)*size(X,2),size(X,3)]));
%% warmstart
% r: the rank of subspace 
% k: the number of Gaussians in MoG
% L: the background of X
% E: the residual of X, E=X-L;
r=2;k=3;
ind=randperm(fix(size(X,2)/4))+50;
X_start=X(:,ind(1:30));
[model]  =warmstart(X_start,r,k);
%% running OMoGMF
model.N=50*size(X,1);model.ro=0.99; 
model.tv.mod=0;model.imgsize=imgsize;model.tv.lamda=1;
tic;
[L,E,F,label,~]= OMoGMF(model,X,3);
toc
%show result
for i=1 :400   
I=[reshape(X(:,i),imgsize),reshape(L(:,i),imgsize),reshape((F(:,i)+0.5),imgsize);
   30*reshape(abs(E(:,i)).*(label(:,i)==3),imgsize),10*reshape(abs(E(:,i)).*(...
   label(:,i)==2),imgsize),2*reshape(abs(E(:,i)).*(label(:,i)==1),imgsize) ];
 imshow(I) ;pause(0.03)
end 
%% running OMoGMF with subsampling
model.lamda=0.99;
model.SSrate=0.01;% 0.01 subsampling rate
tic
[L]=OMoGMF_SS(model,X,3);
toc;
%show result
for i=1 :400   
I=[reshape(X(:,i),imgsize),reshape(L(:,i),imgsize),...
    reshape((X(:,i)-L(:,i))+0.5,imgsize)];
imshow(I) ;pause(0.03)
end 