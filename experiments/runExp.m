cd('..')
addpath(genpath(cd))

[L,S] = generateL(20,5,3,0.1);

%opts.mu = 0.1;
opts.tol = 1e-6;
%opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;
%opts.penalty = 1e4;

[l, s] = testTCPCP(L, S, opts)

l 
s