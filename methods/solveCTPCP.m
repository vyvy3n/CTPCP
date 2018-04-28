addpath(genpath(cd))
clear
% read image
I = double(imread('./image/testimg.jpg'));
X = I/255;
%X = X(21:50,21:50,:);
dim = size(X);
% Noise on X
Xn = X;
%rng(0)
%ind = find(rand(numel(X),1)<0.01);
seed = 0916;
rand('state',seed);
Xn = X+1.5*rand(size(X)).*(rand(size(X))<0.02);

observeNoise = 0; % 1 for "add noise", 0 for "no noise" on observations
lambda = 1/sqrt(max(dim(1:2))*dim(3)) % lambda in "||L||_* + \lambda ||S||_1"

%opts.mu = 0.1;
opts.tol = 1e-6;
%opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;
%opts.penalty = 1e4;

%% Sampling

% q: int, the number of sampled obeservations
q = floor(0.8*numel(Xn));
% Gauss sampling matrix GM
GM = randn(q,numel(Xn));
GM2 = GM'*GM;
g = GM*Xn(:); % sample
% Noise on sampling observations
if observeNoise == 1
    %g = g+ones(q,1).*(0.01*max(abs(g))*normrnd(0,1,[q,1]));
    ind = find(rand(numel(g),1)<0.3);
    g(ind) = rand(length(ind),1);
end
%imshow(reshape(pinv(GM'*GM)*GM'*g,dim))
%% slove TCPCP

tic
[L, S, obj, err, iter] = tcpcp(dim,g,GM,GM2,lambda,opts);
toc

maxP = max(abs(Xn(:)));

figure(1)
subplot(2,2,1)
imshow(X/maxP)% original picture
title('Original Image');
subplot(2,2,2)
imshow(Xn/maxP)% original picture
title('Noise Observation');
subplot(2,2,3)
imshow(L/maxP)% L solved by TCPCP
title('L');
subplot(2,2,4)
imshow(S/maxP)% S solved by TCPCP
title('S');

savefig('TCPCPresult.fig')
save('TCPCPresult.mat')

max(abs(S(:)))