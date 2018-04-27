function [X, L, S] = generateL(n,r,n3,rhoS)

% Generate a rank r (tensor tubal rank) Tensor `L` with shape [n, n, n3]
% and a sparse tensor `S` with the same shape 
% 
% L = P * Q
% 
% where `*` denotes t-product, implemented by 'tensor_tools/tprod.m'
% P: dim = [n, r, n3], i.i.d. sampled from N(0, 1/n) distribution
% Q: dim = [r, n, n3], i.i.d. sampled from N(0, 1/n) distribution
%
% return: L: dim = [n, n, n3], tensor tubal rank of `L` is r
%         S: dim = [n, n ,n3], sampled from i.i.d. Bernoulli +-1 entries

P = randn([n,r,n3])/sqrt(n);
Q = randn([r,n,n3])/sqrt(n);
L = tprod(P, Q);

S = rand([n,n,n3]);
S = 2*(S>0.5)-1; % generate i.i.d. Bernoulli +-1 entries
S = S.*(rand(size(S))<=rhoS); 

% normalization
X = L+S;
maxP = max(abs(X(:)));
X = X/maxP;
L = L/maxP;
S = S/maxP;