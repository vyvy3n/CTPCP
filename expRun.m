addpath(genpath(cd))
clear

Rho = [1.01];
N = [20,10];
Ratio = [0.1];

results = zeros(size(Rho,2)*size(N,2)*size(Ratio,2),11); %placeholders for saving results
i = 0;
for n=N
    for rho=Rho
        for ratio=Ratio
        i = i+1;
        
        r = ratio*n; %tensor nuclear rank = 0.1n, when generate simulation data L
        d = 3;
        rng(1)
        [X,L,S] = generateL(n,r,d,0.02);

        opts.mu = 0.001;
        opts.tol = 1e-6;
        opts.rho = rho;
        opts.max_iter = 1e3*5;
        opts.penalty = 0.5;
        opts.DEBUG = 1;
        opts.muOnly = 1;

        % trpca-tnn
        [L2,S2,obj,err,iter] = trpca_tnn(X,opts);
        errL2 = norm(L(:)-L2(:),2)/norm(L(:),2);
        errS2 = norm(S(:)-S2(:),2)/norm(S(:),2);

        l2 = errL2 < 1e-5;
        s2 = errS2 < 1e-7;

        % tcpcp
        samplingRatio = 0.95; %sampling ratio of TCPCP
        [Lhat, Shat, l, s, errL, errS, runtime] = expTCPCP(L, S, opts, samplingRatio);

        % Plotting
        maxP = 1; 
        figure(i)
        subplot(2,4,1); imshow(L+S/maxP);  title('Origin Image');
        subplot(2,4,2); imshow(L/maxP);    title('L');
        subplot(2,4,6); imshow(S/maxP);    title('S');
        subplot(2,4,3); imshow(Lhat/maxP); title('L_{tcpcp}');
        %subplot(2,4,3); imshow(Lhat/max(max(max(Lhat)))); title('L_{tcpcp}');
        subplot(2,4,7); imshow(Shat/maxP); title('S_{tcpcp}');
        subplot(2,4,4); imshow(L2/maxP);   title('L_{trpca}');
        subplot(2,4,8); imshow(S2/maxP);   title('S_{trpca}');
        % Save results at each step
        expInd = ['n_',num2str(n),'_rho',num2str(rho),'_r',num2str(ratio)];
        saveas(gcf,['TCPCPresult_',expInd,'.png'])
        savefig(['TCPCPresult_',expInd,'.fig'])
        save(['results_',expInd,'.mat'])
        close all

        %l,s(tcpcp),l2,s2(trpac-tnn) indicate if it is a correct recovery, while 1=True
        %errL, errS(tcpcp); errL2, errL2(trpac-tnn) is the relative error
        results(i,1)=l;
        results(i,2)=s;
        results(i,3)=l2;
        results(i,4)=s2;
        results(i,5)=errL;
        results(i,6)=errS;
        results(i,7)=errL2;
        results(i,8)=errS2;
        results(i,9)=n;
        results(i,10)=rho;      
        results(i,11)=runtime;
        results(i,12)=ratio;
        dlmwrite('expResults.csv',results(i, :), 'delimiter',',','-append');
        end
    end
end
save('allExpResults.mat','results')