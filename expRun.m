cd('..')
addpath(genpath(cd))
%rng(0)

[L,S] = generateL(20,20,3,0.1);
opts.mu = 0.01;
opts.tol = 1e-6;
opts.rho = 1.05;
opts.max_iter = 500;
opts.DEBUG = 1;
opts.penalty = 0.0005;

[l, s, errL, errS] = expTCPCP(L, S, opts);
%tcpcp-tnn
X = L+S;
[X2,S2,obj,err,iter] = trpca_tnn(X,opts);
L2 = X2-S2;
errL2 = norm(L(:)-L2(:),2)/norm(L(:),2);
errS2 = norm(S(:)-S2(:),2)/norm(S(:),2);

l2 = errL2 < 1e-5;
s2 = errS2 < 1e-8;

errL
errL2
errS
errS2

maxP = max(abs(X(:)));

figure(1)
subplot(2,3,1)
imshow(L/maxP)
title('L');
subplot(2,3,4)
imshow(S/maxP)
title('S');
subplot(2,3,2)
imshow(l/maxP)%
title('l');
subplot(2,3,5)
imshow(s/maxP)
title('s');
subplot(2,3,3)
imshow(L2/maxP)
title('L2');
subplot(2,3,6)
imshow(S2/maxP)
title('S2');