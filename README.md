# OMoGMF
# OMoGMF and t-OMoGMF
## Robust Online Matrix Factorization for Dynamic Background Subtraction,TPAMI 2017

### Main Contents

**warmstrat**:  some warmstart functions 

**functions**:  some main functions 

**internal**:   some functions for image alignment

**test_data**:     test data:
                   airport.mat   144x176x400  video without camera jitter
                   t_airport.mat 144x176x400  synthetic transformed video with camera jitter
                   sidewalk.mat  220x352x400  real-world video with camera jitter
                   Dummy.mat     59 x 59x100  unaligned face with different illumination

Directly run the testing demo `demo.m` and `demo_t.m` to test OMoGMF and t-OMoGMF,respectively and
run the testing demo `demo_path.m` with your video image path.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the following papers:
% [1] Deyu Meng, Fernando De la Toree, Matrix Factorization with Unknown Noise. ICCV, 2013.
% [2] Hongwei Yong, Deyu Meng, Wangmeng Zuo, Lei Zhang, Robust Online Matrix Factorization for Dynamic Background Subtraction, IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), 2017. In press.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Resources for research purpose only, shall not be used for commercial purposes! All copyrights belong to the original anthors. The technology has applied for patents. If you want to purchase the patents for commercial purposes, please contact the corresponding author: Deyu Meng, dymeng@mail.xjtu.edu.cn. Thank you!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Hongwei Yong. If having any question, feel free to contact: cshyong@comp.polyu.edu.hk or yonghw@stu.xjtu.edu.cn.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.0, release date: 2017.8.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
