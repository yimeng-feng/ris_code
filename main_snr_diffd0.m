%%
% currentFolder = pwd;
% addpath(genpath(currentFolder));
%%
clc, clear
Mt = 4;%BS antenna
Mr =4;%user antenna
Ns=4;%data stream
M =640;%ris element number
p0=16;%ris line number
p1=sqrt(Mr); %user line number
numMC =10;%monte carlo number
snrDbSet =-20:5:10;
c=3e8;%speed of light
fc=7e10;%carrier frequency
lambda=c/fc;%wavelength
d0_1=20*lambda;%distance between ris and BS
d0_2=200*lambda;%distance between ris and BS
d0_3=2000*lambda;%distance between ris and BS
dx_bs=2*lambda;%BS inter-element distance
dx_user=2*lambda;%user inter-element distance
dx_ris=0.5*lambda;%ris inter-element distance
alpha=0.5;%reflection efficiency
Ropt1 = zeros(length(snrDbSet),numMC);
Ropt2= zeros(length(snrDbSet),numMC);
Ropt3= zeros(length(snrDbSet),numMC);
RpsInf1= zeros(length(snrDbSet),numMC);
RpsInf2= zeros(length(snrDbSet),numMC);
RpsInf3= zeros(length(snrDbSet),numMC);
%%
for mm = 1:numMC
    if mod(mm,10)==1
        mm
    end
    mm
    % generate chanel matrix
    %channel between RIS to user
    Nc=3;%cluster number
    Nray=5;%ray number
    G = channel_generation_ura(M,Mr,Nc,Nray,dx_ris,dx_user,lambda,p0,p1);
    H1=channel_generation(M,Mt,d0_1,alpha,p0,lambda,dx_bs,dx_ris);
    H2=channel_generation(M,Mt,d0_2,alpha,p0,lambda,dx_bs,dx_ris);
    H3=channel_generation(M,Mt,d0_3,alpha,p0,lambda,dx_bs,dx_ris);
    Hd= zeros(Mr,Mt);%the direction channel between BS and User
    %% calculate rates
    for indxSnrDb = 1:length(snrDbSet)
        snrDb = snrDbSet(indxSnrDb);
        snrLin = db2pow(snrDb);
        %% upper bound 20lambda
        Ropt1(indxSnrDb,mm) = RIS_upperbound(G,H1,snrLin,Ns);
        %% perfect PS based, No Switch 20lambda
        RpsInf1(indxSnrDb,mm) = RIS(G,H1,snrLin,Ns,Hd,'psInf');
        %% upper bound 200lambda
        Ropt2(indxSnrDb,mm) = RIS_upperbound(G,H2,snrLin,Ns);
        %% perfect PS based, No Switch 200lambda
        RpsInf2(indxSnrDb,mm) = RIS(G,H2,snrLin,Ns,Hd,'psInf');
        %% upper bound 2000lambda
        Ropt3(indxSnrDb,mm) = RIS_upperbound(G,H3,snrLin,Ns);
        %% perfect PS based, No Switch 2000lambda
        RpsInf3(indxSnrDb,mm) = RIS(G,H3,snrLin,Ns,Hd,'psInf');
    end
end
%% det(Popt'*chanMat'*chanMat*Popt)
figure
width = 1.5;

plot(snrDbSet,mean(RpsInf1,2),'k-d','LineWidth',width), hold on
plot(snrDbSet,mean(Ropt1,2),'k','LineWidth',width), hold on
plot(snrDbSet,mean(RpsInf2,2),'b--d','LineWidth',width), hold on
plot(snrDbSet,mean(Ropt2,2),'b','LineWidth',width), hold on
plot(snrDbSet,mean(RpsInf3,2),'r-d','LineWidth',width), hold on
plot(snrDbSet,mean(Ropt3,2),'r','LineWidth',width), hold on
legend('Our proposed algorithm (d_0=20\lambda)','The Upper Bound (d_0=20\lambda)','Our proposed algorithm (d_0=100\lambda)','The Upper Bound (d_0=100\lambda)','Our proposed algorithm (d_0=1000\lambda)','The Upper Bound (d_0=1000\lambda)','Location','NorthWest')%'RIS scheme (b=1)','RIS scheme (b=2)',
grid on
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bps/Hz)')
title(['URA in a Channel with ' num2str(Nc*Nray) ' Multipath Clusters, M_t=' num2str(Mt) ', M_r=' num2str(Mr) ', N_s=' num2str(Ns) ', M=' num2str(M) ', \beta=' num2str(alpha)]);
return
