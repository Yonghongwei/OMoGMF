% This is the testing demo for video background subtraction with OMoGMF on your video path

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
para.k=3;
para.r=2;
para.display=1;
para.ro=0.98;
para.N=50;
para.startindex=[200,1000];
para.startnumber=30;
para.iter=3;
video_path='G:\BS Videos\Li\Ariport\';%add your video image path 
run_video(video_path,para)