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
%tcpcp-tnn
[X2,S2,obj,err,iter] = trpca_tnn(Xn,lambda,opts);
L2 = X2-S2;
errL2 = norm(L(:)-L2(:),2)/norm(L(:),2);
errS2 = norm(S(:)-S2(:),2)/norm(S(:),2);

l2 = errL2 < 1e-5;
s2 = errS2 < 1e-8;

errL
errL2
errS
errS2