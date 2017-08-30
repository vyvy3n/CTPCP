addpath(genpath(cd))
clear

pic_name = [ './image/testimg.jpg'];
I = double(imread(pic_name));
X = I/255;
    
[n1,n2,n3] = size(X);

opts.mu = 1e-3;
opts.tol = 1e-6;
opts.rho = 1.2;
opts.max_iter = 500;
opts.DEBUG = 1;

% p = 0.6;
% maxP = max(abs(X(:)));
% omega = find(rand(n1*n2*n3,1)<p);
% M = zeros(n1,n2,n3);
% M(omega) = X(omega);

%% Sampling 'matrix'(tensors)
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
%% 
