% This is the testing demo for video background subtraction and face alignment with t-OMoGMF and it-OMoGMF
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

data={'test_data\t_airport.mat','test_data\sidewalk.mat','test_data\Dummy.mat'};

for data_ind=1:3

if data_ind==1
%% synthetic  data: airport
load(data{data_ind});
imgsize=[size(X,1),size(X,2)];
X=im2double(reshape(X,[size(X,1)*size(X,2),size(X,3)]));

%warmstart
ind=randperm(fix(size(X,2)/3))+40;
X_start=X(:,ind(1:30));
X_start=reshape(X_start,[imgsize,size(X_start,2)]);
[model] =t_warmstart(X_start,2,2);

%online running
model.N=50*size(X,1);model.ro=0.99;
model.imgsize=imgsize;model.preAlignment=1;
Tau=zeros(6,size(X,2));iter=10;
[L,B,E,label,model,Tau]=  t_OMoGMF(model,X(:,1:400),iter,Tau);

%show result
for i=1:400
imshow([reshape(X(:,i),imgsize),reshape(B(:,i),imgsize);
    reshape(L(:,i),imgsize),reshape(abs(E(:,i)),imgsize)]);  
pause(0.01)   
end

elseif data_ind==2
%% real-world data: sidewalk
load(data{data_ind});
imgsize=[size(X,1),size(X,2)];
X=im2double(reshape(X,[size(X,1)*size(X,2),size(X,3)]));
%warmstart
ind=randperm(fix(size(X,2)/3));
X_start=X(:,ind(1:30));
X_start=reshape(X_start,[imgsize,size(X_start,2)]);
[model] =t_warmstart(X_start,2,2);
%online running
model.N=50*size(X,1);model.ro=0.99;
model.imgsize=imgsize;model.preAlignment=1;%if not use prealignment, set a large iter
Tau=zeros(6,size(X,2));iter=5;
[L,B,E,label,model,Tau]= t_OMoGMF(model,X(:,1:400),iter,Tau);
%show result
for i=1:400
imshow([reshape(X(:,i),imgsize),reshape(B(:,i),imgsize);
    reshape(L(:,i),imgsize),reshape(abs(E(:,i)),imgsize)]);  
pause(0.01)   
end
elseif data_ind==3
%% face data: Dummy
load(data{data_ind});
imgsize=[size(X,1),size(X,2)];
X=im2double(reshape(X,[size(X,1)*size(X,2),size(X,3)]));

%warmstart
ind=randperm(fix(size(X,2)));
X_start=X(:,ind(1:50));
X_start=reshape(X_start,[imgsize,size(X_start,2)]);
[model] =t_warmstart(X_start,5,2);
%batch running
model.N=50*size(X,1);model.ro=0.98;
model.imgsize=imgsize;model.preAlignment=0;%if not use prealignment, set a large iter
iter=5;
[model,label,L,E,TX,Tau]=it_OMoGMF(X,model,iter);
%show result
I0{1}=X;I0{2}=TX;I0{3}=L;I0{4}=abs(E.*(label==1))*2;I0{5}=abs(E.*(label==2))*10;
titles={'Original images','Alinged images','Low rank component','Large variance residual','Small variance residual'};
for k=1:5
    I=reshape(I0{k},[imgsize,10,10]);
for i=1:10;for j=1:10 Im(((i-1)*size(I,1)+1):((i)*size(I,1))...
         ,((j-1)*size(I,2)+1):((j)*size(I,2)))=I(:,:,i,j); end; end
figure(k)
imshow(Im)
title(titles{k})
end
    
end
end
 
 