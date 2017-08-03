function  run_video(video_path,para)
%%online mog matrix factorization for video background subtraction
%input: video_path the video images path
%       para  all  hyperparameters for OMoGMF model 
%       para.N controls the speed of updating MoG parameters
%       para.ro controls the speed of updating U
%       para.tv.mod  switch of TV threshold: if tv.mod==1 using TV threshold
%       para.tv.lamda  parameters of TV threshold
%       para.iter the number of iterations
%       para.r the rank of subspace
%       para.k the number of Gaussians in MoG
%       para.display  if display==1, show the result
%       para.startindex the  frame index range [startindex(1),startindex(2)] for warmstart
%       para.startnumber  the number of frames for warmstart
%Written by Hongwei Yong(cshyong@comp.polyu.edu.hk or yonghw@stu.xjtu.edu.cn).
%% init
path=dir(video_path);
len=length(path);
if  nargin<2
para=init(len);
else
para=init(len,para);   
end
k=para.k;
r=para.r;
iter=para.iter;
display=para.display;
%% warmstart 
startindex=para.startindex;
startnumber=para.startnumber;

im=imread([video_path,path(3).name]);
imgsize=[size(im,1),size(im,2),size(im,3)];
imgnum_start=startindex(2)-startindex(1)+1;
ind=randperm(imgnum_start);

X_start=zeros(imgsize(1)*imgsize(2)*imgsize(3),startnumber);
for i=1:startnumber
im=im2double(imread([video_path,path(ind(i)+2).name]));
X_start(:,i)=im(:);
end

[model]  =warmstart(X_start,r,k);
clear X_start
model.ro=para.ro;
model.N=para.N*imgsize(1)*imgsize(2)*imgsize(3);
model.tv.mod=para.tv.mod;
model.tv.lamda=para.tv.lamda;
%% main 

for i=1:len-2
   if mod(i,100)==0||i==1
      disp(['Calculating the model of the ',num2str(i),'th frame']);
   end
   im=im2double(imread([video_path,path(i+2).name]));
   x=im(:);
  [model,v,label,~] =onlinemogmf(model,x,iter,ones(size(x)));% main function
   U=[model.U];
   L=U*v;E=x-L;
  F0=E.*(label==1);
 % TV threshood
 if model.tv.mod==1
  FF=reshape(E,imgsize);
  lam = model.tv.lamda*(model.Sigma(1));
  F0= TVthre(FF,lam);
 end
 %% show
F=F0(:);
if display==1
I=[reshape(x,imgsize),reshape(L,imgsize),reshape((E+0.5),imgsize);
   30*reshape(abs(E).*(label==3),imgsize),10*reshape(abs(E).*(...
   label==2),imgsize),2*reshape(abs(E).*(label==1),imgsize) ];
 imshow(I) ;pause(0.001)
end
 end
%% output
end
 
function para=init(len,para)
if nargin<2
  para.r=3;para.k=3;
end
if (~isfield(para,'r'))
para.r=2;
end
if (~isfield(para,'k'))
para.k=3;
end
if (~isfield(para,'display'))
para.display=1;
end
if (~isfield(para,'ro'))
para.ro=0.985;
end
if (~isfield(para,'N'))
para.N=50;
end
if (~isfield(para,'iter'))
para.iter=3;
end

if (~isfield(para,'startindex'))
para.startindex=[1,len-2];
end

if (~isfield(para,'startnumber'))
para.startnumber=50;
end

if (~isfield(para,'tv'))
       para.tv.mod=0;
       para.tv.lamda=1;
end
if (~isfield(para.tv,'.mod'))
   para.tv.mod=0;
end
if (~isfield(para.tv,'lamda'))
   para.tv.lamda=1;
end
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
