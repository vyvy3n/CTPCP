addpath(genpath(cd))
clear

I = double(imread('./image/testimg.jpg'));
X = I/255;
X = X(1:50,1:50,:);
dim = size(X);

[n1,n2,n3] = size(X);
noise = 0; % 1 for "add noise", 0 for "no noise" on observations
lambda = 1;

opts.mu = 1e-3;
opts.tol = 1e-6;
opts.rho = 1.2;
opts.max_iter = 500;
opts.DEBUG = 1;

%% Sampling
% q: int, the number of sampled obeservations
q = floor(0.5*numel(X));
% % Gauss Sampling Tensors
% G = cell(q,1);
% % sampling
% g = zeros(q,1);
% for i=1:q
%     G{i} = normrnd(0,1,size(X));
%     g(i) = dot(G{i}(:),X(:));
% end

% Gauss sampling matrix GM
GM = randn(q,numel(X));
g = GM*X(:);
%% Noise
if noise == 1
    g = g+ones(q,1).*(0.01*max(abs(g))*normrnd(0,1,[q,1]));
end

%% slove TCPCP

[Xhat,err,iter] = tcpcp(dim,g,GM,lambda,opts);

err
iter

maxP = max(abs(X(:)));

figure(1)
subplot(1,2,1)
imshow(X/maxP)
% subplot(1,3,2)
% imshow(M/maxP)
subplot(1,2,2)
imshow(Xhat/maxP)