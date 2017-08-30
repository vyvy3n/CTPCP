addpath(genpath(cd))
clear

pic_name = [ './image/testimg.jpg'];
I = double(imread(pic_name));
X = I/255;
% X = X(1:10,1:10,3);

[n1,n2,n3] = size(X);
noise = 1 % 1 for "add noise", 0 for "no noise" on observations

opts.mu = 1e-3;
opts.tol = 1e-6;
opts.rho = 1.2;
opts.max_iter = 500;
opts.DEBUG = 1;

%% Sampling
% q: int, the number of sampled obeservations
q = floor(0.1*numel(X));
% Gauss Sampling Tensors
G = cell(q,1);
% sampling
GM = zeros(q,1);
for i=1:q
    G{i} = normrnd(0,1,size(X));
    GM(i) = dot(G{i}(:),X(:));
end

%% Noise
if noise == 1
    GM = GM+ones(q,1).*(0.01*max(abs(GM))*normrnd(0,1,[q,1]));
end

%% slove TCPCP

[Xhat,err,iter] = tcpcp(GM,opts);

err
iter

figure(1)
subplot(1,2,1)
imshow(X/maxP)
% subplot(1,3,2)
% imshow(M/maxP)
subplot(1,2,2)
imshow(Xhat/maxP)
