cd('..')
addpath(genpath(cd))

[L,S] = generateL(20,15,3,0.1);
opts.mu = 0.01;
opts.tol = 1e-6;
opts.rho = 1.05;
opts.max_iter = 500;
opts.DEBUG = 1;
opts.penalty = 0.0005;

[l, s, errL, errS] = expTCPCP(L, S, opts);

l 
s