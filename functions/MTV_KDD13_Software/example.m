%% this is an example
% parameter setting
addpath('Code');
DISPLAY_SETTING = 0; % display the primal & dual error 
MAX_ITER = 100;     % the maximum iteration number
TOL_SETTING = [1e-3, 1e-3]; % the absolute tolerance and relative tolerance
lam = 0.35;          % regularization paramter
rho = 10;            % the dual update length for ADMM

N = 200;             % the size of image N x N

I=MakeRDSquares(N,7,40); % generate the image with 7 blocks
Im = I + 0.2*randn(N,N);
% set initial solution
[X,iter] = mtv2(Im,lam,rho,MAX_ITER,TOL_SETTING,DISPLAY_SETTING); %iter is number of iteration

[Xp,iterp] = mtvp2(Im,lam,rho,MAX_ITER,TOL_SETTING,DISPLAY_SETTING); % parallel version
figure(1);
subplot(1,3,1);
imshow(I);title('Ground Truth');
subplot(1,3,2);
imshow(Im);title('Noisy image');
subplot(1,3,3);
imshow(X);title('Recovered');

% generate 3D images
I3 = Make3DImage(N,51);
Im3 = I3 + 0.2*randn(N,N,51);

[X3,iterp3]=mtvp3(Im3, lam,rho,MAX_ITER,TOL_SETTING,DISPLAY_SETTING);

figure(2);
subplot(1,3,1);
imshow(I3(:,:,1));title('Ground Truth');
subplot(1,3,2);
imshow(Im3(:,:,1));title('Noisy image');
subplot(1,3,3);
imshow(X3(:,:,1));title('Recovered');
