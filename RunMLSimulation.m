%% Estimation of bias in return values, exceedance probabilities and log exceedance probabilities 
% Generalised Pareto tail fitting
% ML, MOM, PWM, EB methods

close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
clc; clear;

%% STEP 1 : Generate large set of parameter estimates

Xi=(-0.4:0.05:0.1)';
N=10.^(2:4)';
nXi=size(Xi,1);
nN=size(N,1);
nRpl=100; %Set this to at least 1e4. Set at 100 so that a new user can see full analysis quickly for any bugs
Est=1; %estimation method

if 1; %Calculate true return values, and generate sample estimates
    
    %Generate basic estimates for xi and sgm
    fprintf(1,'Total number of fits is %g\n',nXi*nRpl*nN);
    iC=0;
    P=nan(nN,nXi,nRpl,2);
    for iN=1:nN;
        for iX=1:nXi;
            for iR=1:nRpl;
                iC=iC+1;
                X=gprnd(Xi(iX)*ones(N(iN),1),1,0);
                switch Est;
                    case 1; %ML
                        P(iN,iX,iR,:)=gpfit(X);
                    case 2; %MOM
                        P(iN,iX,iR,:)=gpfitMOM(X);
                    case 3; %PWM
                        P(iN,iX,iR,:)=gpfitPWM(X);
                    case 4; %EB (Zhang)
                        P(iN,iX,iR,:)=gpfitEB(X);
                end;
                if rem(iC,100)==0; fprintf(1,'%g\n',iC); end;
            end;
        end;
    end;
    
    tStr=sprintf('Data'); %given different name to avoid overwrite
    save(tStr,'P','Xi','N','nXi','nN','nRpl');
    
else;
    
    tStr=sprintf('Datam0p4t0p1'); %given different name to avoid overwrite
    %tStr=sprintf('DatGold'); %given different name to avoid overwrite
    load(tStr);
    
end;

%% STEP 2 : Specify Rat and RP and calculate bias in return values, exceedance probabilities and log exceedance probabilities using Maximum Likelihood

Rat=100; %annual expected rate of occurrence
RP=10000; %return period

if 1; %Estimate return values, exceedance probabilities and log exceedance probabilities using Maximum Likelihood
    
    %Calculate true return values
    for iN=1:nN;
        for iX=1:nXi;
            E0(iN,iX)=RtrVal(Xi(iX),1,Rat*RP);
            T0(iN,iX)=AnnExcPrb(E0(iN,iX),Xi(iX),1,Rat);
        end;
    end;
    F0=E0.^2;
    
    E=nan(nN,nXi,4);
    T=nan(nN,nXi,4);
    F=nan(nN,nXi,4);
    G=nan(nN,nXi,4);
    if 1; %Estimate return values and probabilities
        
        EXi=permute(mean(P,3),[1 2 4 3]);
        for iN=1:nN;
            for iX=1:nXi;
                
                fprintf(1,'N=%g,Xi=%g\n',N(iN),Xi(iX));
                
                % E1=Q_A(1-1/N|E[P])
                E(iN,iX,1)=RtrVal(EXi(iN,iX,1),EXi(iN,iX,2),Rat*RP);
                T(iN,iX,1)=AnnExcPrb(E(iN,iX,1),Xi(iX),1,Rat);
                F(iN,iX,1)=RtrVal(EXi(iN,iX,1),EXi(iN,iX,2),Rat*RP).^2;
                G(iN,iX,1)=AnnExcPrb(F(iN,iX,1).^0.5,Xi(iX),1,Rat);
                
                % E2=E[Q_A(1-1/N|P)]
                t1=nan(nRpl,1);
                t2=nan(nRpl,1);
                for iR=1:nRpl;
                    t1(iR)=RtrVal(P(iN,iX,iR,1),P(iN,iX,iR,2),Rat*RP);
                    t2(iR)=RtrVal(P(iN,iX,iR,1),P(iN,iX,iR,2),Rat*RP).^2;
                end;
                E(iN,iX,2)=mean(t1);
                T(iN,iX,2)=AnnExcPrb(E(iN,iX,2),Xi(iX),1,Rat);
                F(iN,iX,2)=mean(t2);
                G(iN,iX,2)=AnnExcPrb(F(iN,iX,2).^0.5,Xi(iX),1,Rat);
                
                %E3 = QPP_A_N(exp(-1))
                DInc=100;DMxm=10000;
                SlnLwr=0;
                for iI=1:4;
                    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
                    tCdf=AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,RP)';
                    SlnLwr=D(tCdf<exp(-1));SlnLwr=SlnLwr(end);
                    DInc=DInc/100; DMxm=DMxm/100;
                end;
                E(iN,iX,3)=SlnLwr;
                T(iN,iX,3)=AnnExcPrb(E(iN,iX,3),Xi(iX),1,Rat);
                %
                DInc=100;DMxm=10000;
                SlnLwr=0;
                for iI=1:4;
                    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
                    tCdf=AnnCdf2(D.^0.5,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,RP)';
                    SlnLwr=D(tCdf<exp(-1));SlnLwr=SlnLwr(end);
                    DInc=DInc/100; DMxm=DMxm/100;
                end;
                F(iN,iX,3)=SlnLwr;
                G(iN,iX,3)=AnnExcPrb(F(iN,iX,3).^0.5,Xi(iX),1,Rat);
                
                %E4 = QPP_A(1-1/N)
                DInc=100;DMxm=10000;
                SlnLwr=0;
                for iI=1:4;
                    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
                    tCdf=AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,1)';
                    SlnLwr=D(tCdf<1-1/RP);SlnLwr=SlnLwr(end);
                    DInc=DInc/100; DMxm=DMxm/100;
                end;
                E(iN,iX,4)=SlnLwr;
                T(iN,iX,4)=AnnExcPrb(E(iN,iX,4),Xi(iX),1,Rat);
                %
                DInc=100;DMxm=10000;
                SlnLwr=0;
                for iI=1:4;
                    D=(0:DInc:DMxm-DInc)'+SlnLwr; nD=size(D,1);
                    tCdf=AnnCdf2(D.^0.5,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,1)';
                    SlnLwr=D(tCdf<1-1/RP);SlnLwr=SlnLwr(end);
                    DInc=DInc/100; DMxm=DMxm/100;
                end;
                F(iN,iX,4)=SlnLwr;
                G(iN,iX,4)=AnnExcPrb(F(iN,iX,4).^0.5,Xi(iX),1,Rat);
                
                %E5 0.75E3+0.25E4
                E(iN,iX,5)=0.5*E(iN,iX,3)+0.5*E(iN,iX,4);
                T(iN,iX,5)=AnnExcPrb(E(iN,iX,5),Xi(iX),1,Rat);
            end;
        end;
    end;
    
    BE=(E-E0)./E0;
    BT=(T-T0);
    BG=G-T0; %T0 and G0 must be the same
    BF=(F-F0)./F0;
    
    tStr=sprintf('BiasRat%gRP%g',Rat,RP);
    save(tStr,'E','T','E0','T0','BE','BT','Xi','nXi','N','nN','Rat','RP','F0','F','BF','G','BG');
    
else; %load previous work
    
    tStr=sprintf('BiasRat%gRP%g',Rat,RP);
    load(tStr,'E','T','E0','T0','BE','BT','Xi','nXi','N','nN','Rat','RP','F0','F','BF','G','BG');
    
end;

%% STEP 3 : Standard plots of results

if 1; 
    figure(1); clf;
    ty=BE(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Fractional bias of return value');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$'}; pLgn(tLgn,'NW',30,3);
    tStr=sprintf('FractionalBiasReturnValue%g',RP);
    pGI(tStr,2);
    
    figure(2); clf;
    ty=BT(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of exceedance probability of return value');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$'}; pLgn(tLgn,'NE',30,3);
    tStr=sprintf('BiasExcProb%g',RP);
    pGI(tStr,2);
    
    figure(3); clf;
    ty=E(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,E0(iN,:)','color',ones(3,1)*0.9,'linewidth',15);
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Return value');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$'}; pLgn(tLgn,'NW',30,3);
    tStr=sprintf('ReturnValue%g',RP);
    pGI(tStr,2);
    
    figure(4); clf;
    ty=BF(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Fractional bias of return value');
    tLgn={'$\hat{q}_1^*$';'$\hat{q}_2^*$';'$\hat{q}_3^*$';'$\hat{q}_4^*$'}; pLgn(tLgn,'NW',30,3);
    tStr=sprintf('FractionalBiasReturnValueForResponse%g',RP);
    pGI(tStr,2);
    
    figure(5); clf;
    ty=F(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,F0(iN,:)','color',ones(3,1)*0.9,'linewidth',15);
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Return value of response');
    tLgn={'$\hat{q}_1^*$';'$\hat{q}_2^*$';'$\hat{q}_3^*$';'$\hat{q}_4^*$'}; pLgn(tLgn,'NW',30,3);
    tStr=sprintf('ReturnValueForResponse%g',RP);
    pGI(tStr,2);
    
    if exist('BG')==1;
        figure(6); clf;
        ty=BG(:,:,[1 2 4 3]);
        for iN=1:nN;
            subplot(1,3,iN); hold on;
            for j=1:4;
                plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
            end;
            tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
            pDflBig;pAxsLmt;
        end;
        subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of exceedance probability for response');
        tLgn={'$\hat{q}_1^*$';'$\hat{q}_2^*$';'$\hat{q}_3^*$';'$\hat{q}_4^*$'}; pLgn(tLgn,'NE',30,3);
        tStr=sprintf('BiasExcProbForResponse%g',RP);
        pGI(tStr,2);
    end;
    
    figure(7); clf;
    t=log10(T./T0); t(isinf(t))=nan;
    ty=t(:,:,[1 2 4 3]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of $\log_{10}$ exceedance probability','Interpreter','latex');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$'}; pLgn(tLgn,'SW',30,3);
    tStr=sprintf('BiasLogExcProb%g',RP);
    pGI(tStr,2);
    
    figure(8); clf;
    ty=BE(:,:,[1 2 4 3 5]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,ty(iN,:,5)','linestyle',pLinStl(1),'color',ones(3,1)*0.8,'linewidth',10,'marker','none');
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Fractional bias of return value');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$';'$(\hat{q}_3+\hat{q}_4)/2$'};
    pLgn(tLgn,'NW',30,3,...
        cellstr(char(pLinStl(1),pLinStl(2),pLinStl(3),pLinStl(4),pLinStl(1))),...
        {ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0.8;});
    tStr=sprintf('FractionalBiasReturnValuePlus%g',RP);
    pGI(tStr,2);
    
    figure(9); clf;
    ty=BT(:,:,[1 2 4 3 5]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,ty(iN,:,5)','linestyle',pLinStl(1),'color',ones(3,1)*0.8,'linewidth',10,'marker','none');
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of exceedance probability');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$';'$(\hat{q}_3+\hat{q}_4)/2$'};
    pLgn(tLgn,'NE',30,3,...
        cellstr(char(pLinStl(1),pLinStl(2),pLinStl(3),pLinStl(4),pLinStl(1))),...
        {ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0.8;});
    tStr=sprintf('BiasExcProbPlus%g',RP);
    pGI(tStr,2);
    
    figure(10); clf;
    t=log10(T./T0); t(isinf(t))=nan;
    ty=t(:,:,[1 2 4 3 5]);
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,ty(iN,:,5)','linestyle',pLinStl(1),'color',ones(3,1)*0.8,'linewidth',10,'marker','none');
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of $\log_{10}$ exceedance probability','Interpreter','latex');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$';'$(\hat{q}_3+\hat{q}_4)/2$'};
    pLgn(tLgn,'MW',30,3,...
        cellstr(char(pLinStl(1),pLinStl(2),pLinStl(3),pLinStl(4),pLinStl(1))),...
        {ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0.8;});tStr=sprintf('BiasLogExcProbPlus%g',RP);
    pGI(tStr,2);
    
end;

%% STEP 4 : Other diagnostic plots

if 0; % RMSE estimates
    
    figure(1); clf;
    
    for iN=1:nN;
        subplot(1,3,iN); hold on;
        plot(Xi,ty(iN,:,5)','linestyle',pLinStl(1),'color',ones(3,1)*0.8,'linewidth',10,'marker','none');
        for j=1:4;
            plot(Xi,ty(iN,:,j)','linestyle',pLinStl(j),'color','k','linewidth',3,'marker','o');
        end;
        tStr=sprintf('Sample size = %g\n',N(iN)); title(tStr);
        pDflBig;pAxsLmt;
    end;
    subplot(1,3,1); hold on; xlabel('$\xi_0$','interpreter','latex'); ylabel('Bias of $\log_{10}$ exceedance probability','Interpreter','latex');
    tLgn={'$\hat{q}_1$';'$\hat{q}_2$';'$\hat{q}_3$';'$\hat{q}_4$';'$(\hat{q}_3+\hat{q}_4)/2$'};
    pLgn(tLgn,'MW',30,3,...
        cellstr(char(pLinStl(1),pLinStl(2),pLinStl(3),pLinStl(4),pLinStl(1))),...
        {ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0;ones(3,1)*0.8;});tStr=sprintf('BiasLogExcProbPlus%g',RP);
    pGI(tStr,2);
    
end;

if 0; % check inequalities
    Inq=nan(nN,nXi,3);
    for iN=1:nN;
        for iX=1:nXi;
            t=permute(P(iN,iX,:,:),[3 4 1 2]);
            if sum(t(:,1)>0)==0;
                Inq(iN,iX,1)=sum(t(:,2))/(-sum(t(:,1)))-1/(-Xi(iX));
                Inq(iN,iX,2)=sum(t(:,2)./(-t(:,1)))/size(t,1)-1/(-Xi(iX));
                Inq(iN,iX,3)=max(t(:,2)./(-t(:,1)))-1/(-Xi(iX));
            end;
        end;
    end;
end;

if 0; % parameter estimate scatter plots
    clf;
    iX=5;
    for j=1:3;
        subplot(2,3,j);
        iN=j;
        plot(permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),'k.','markersize',10);
        pAxsLmt; pDflBig;
        tStr=sprintf('n=%g, $\\xi_0$=%g',N(iN),Xi(iX));
        title(tStr);
    end;
    iX=11;
    for j=1:3;
        subplot(2,3,3+j);
        iN=j;
        plot(permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),'k.','markersize',10);
        pAxsLmt; pDflBig;
        tStr=sprintf('n=%g, $\\xi_0$=%g',N(iN),Xi(iX));
        title(tStr);
    end;
    subplot(2,3,4); xlabel('$\hat{\xi}$','interpreter','latex'); ylabel('$\hat{\sigma}$');
    pGI('SctPrmEst',2);
end;


if 0; % comparison of predictive distributions
    RP=100;
    clf;
    subplot(1,2,1); hold on;
    iN=2; iX=5;
    DInc=0.01;
    DMnm=1.5;
    DMxm=5.5;
    D=(DMnm:DInc:DMxm)'; nD=size(D,1);
    t0=1-AnnExcPrb(D,Xi(iX),1,Rat);
    t0N=1-AnnExcPrb(D,Xi(iX),1,Rat*RP);
    t1=AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,1)';
    t2=(AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,RP)');
    t3=t2.^(1/RP);
    plot(D,t0,'color',ones(3,1)*0.5,'linestyle',pLinStl(1),'linewidth',3);
    plot(D,t0N,'color',ones(3,1)*0.8,'linestyle',pLinStl(1),'linewidth',3);
    plot(D,t1,'k','linestyle',pLinStl(3),'linewidth',3);
    plot(D,t2,'k','linestyle',pLinStl(4),'linewidth',3,'color',ones(3,1)*0.7);
    plot(D,t3,'k','linestyle',pLinStl(4),'linewidth',3);
    pAxsLmt; pDflBig;
    pLgn({'$F_A$';'$F_{A_N}$';'$\tilde{F}_A$';'$\tilde{F}_{A_N}$';'$\tilde{F}_{A_N}^{1/N}$'},'NW',28,3,...
        cellstr(char(pLinStl(1),pLinStl(1),pLinStl(3),pLinStl(4),pLinStl(4))),...
        {ones(3,1)*0.5;ones(3,1)*0.8;zeros(3,1);ones(3,1)*0.7;zeros(3,1);});
    xlabel('$x$','interpreter','latex');
    ylabel('Cumulative distribution function (cdf)');
    
    subplot(2,2,2); cla; hold on;
    DInc=0.01;
    DMnm=1.5;
    DMxm=5.5;
    D=(DMnm:DInc:DMxm)'; nD=size(D,1);
    t0=1-AnnExcPrb(D,Xi(iX),1,Rat);
    t0N=1-AnnExcPrb(D,Xi(iX),1,Rat*RP);
    t1=(AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,1));
    t2=(AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,RP));
    t3=t2.^(1/RP);
    IsOK=D>=0;
    plot(D(IsOK),(t1(IsOK)-t0(IsOK))*5,'color',ones(3,1)*0,'linewidth',3,'marker','none');
    IsOK=D>=0;
    plot(D(IsOK),t2(IsOK)-t0N(IsOK),'color',ones(3,1)*0.8,'linewidth',3,'marker','none');
    IsOK=D>=0;
    plot(D(IsOK),(t3(IsOK)-t1(IsOK))/5,'color',ones(3,1)*0.5,'linewidth',3,'marker','none');
    pAxsLmt; pDflBig;
    pLgn({'$(\tilde{F}_{A}-F_{A})\times5$';'$\tilde{F}_{A_N}-F_{A_N}$';},'NW',28,3,...
        cellstr(char(pLinStl(1),pLinStl(1))),...
        {ones(3,1)*0;ones(3,1)*0.8;},...
        0.12);
    %Custom legend
    limits = objbounds(findall(gca));
    xl=limits(1); xh=limits(2); yl=limits(3); yh=limits(4); zl=limits(5); zh=limits(6);
    rx=xh-xl; ry=yh-yl; j=1; Gap=0.12; nT=1;
    plot([xl xl+0.2*rx],[yl yl]+Gap*(nT-j+1)*ry,'linestyle','-','linewidth',3,'marker','o','color',ones(3,1)*0.5);
    text(xl+0.25*rx,yl+Gap*(nT-j+1)*ry,'$(\tilde{F}_{A_N}^{1/N}-\tilde{F}_{A_N})/5$','FontSize',28,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
    %True return value
    t=get(gca,'ylim');
    q_0=(1./Xi(iX))*((Rat*RP)^Xi(iX)-1);
    plot(q_0*ones(2,1),t,'k-');
    ylabel('Difference in cdf');
    
    subplot(2,2,4); cla; hold on;
    DInc=0.01;
    DMnm=3.9;
    DMxm=4.5;
    D=(DMnm:DInc:DMxm)'; nD=size(D,1);
    t0=1-AnnExcPrb(D,Xi(iX),1,Rat);
    t0test=AnnCdf2(D,Xi(iX)*ones(10,1),ones(10,1),Rat,1)';
    t1=AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,1);
    t2=(AnnCdf2(D,permute(P(iN,iX,:,1),[3 1 2 4]),permute(P(iN,iX,:,2),[3 1 2 4]),Rat,RP)).^(1/RP);
    IsOK=D>=3.9 & D<=4.5 & t0>=0.98;
    plot(D(IsOK),t0(IsOK),'color',ones(3,1)*0.7,'linestyle',pLinStl(1),'linewidth',3);
    IsOK=D>=3.9 & D<=4.5 & t1>=0.98;
    plot(D(IsOK),t1(IsOK),'k','linestyle',pLinStl(3),'linewidth',3);
    IsOK=D>=3.9 & D<=4.5 & t2>=0.98;
    plot(D(IsOK),t2(IsOK),'k','linestyle',pLinStl(4),'linewidth',3);
    pAxsLmt; pDflBig;
    pLgn({'$F_A$';'$\tilde{F}_A$';'$\tilde{F}_{A_N}^{1/N}$'},'NW',28,3,...
        cellstr(char(pLinStl(1),pLinStl(3),pLinStl(4))),...
        {ones(3,1)*0.5;zeros(3,1);zeros(3,1)},...
        0.12)
    t=get(gca,'xlim')';
    plot(t,(1-1/RP)*ones(2,1),'k-');
    t=get(gca,'ylim');set(gca,'ylim',[t(1) 1]);
    ylabel('cdf');
    
    pGI('q3a4Eda',2);
end;