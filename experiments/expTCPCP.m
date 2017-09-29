function [l, s, errL, errS] = expTCPCP(L, S, opts)

X = L+S;
dim = size(X);
lambda = 1/sqrt(max(dim(1:2))*dim(3)); % lambda in "||L||_* + \lambda ||S||_1"

%% Sampling

% q: int, the number of sampled obeservations
q = floor(0.9*numel(X));
% Gauss sampling matrix GM
GM = randn(q,numel(X));
GM2 = GM'*GM;
g = GM*X(:); % sample

%% slove TCPCP

tic
[Lhat, Shat, obj, err, iter] = tcpcp(dim,g,GM,GM2,lambda,opts);
toc

errL = norm(L(:)-Lhat(:),2)/norm(L(:),2);
errS = norm(S(:)-Shat(:),2)/norm(S(:),2);
l = errL < 1e-5;
s = errS < 1e-8;
