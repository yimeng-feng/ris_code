%%
% currentFolder = pwd;
% addpath(genpath(currentFolder));
%%
clc, clear
Mt = 4;%BS antenna
Mr = 4;%user antenna
Ns = 4;%data stream
M = 640;%ris element number
p0 = 16;%ris line number
p1 = sqrt(Mr); %user line number
numMC =100;%monte carlo number
snrDbSet =-20:5:10;
chanType = 'URA';
c=3e8;%speed of light
fc=7e10;%carrier frequency
lambda=c/fc;%wavelength
d0=10000*lambda;%distance between ris and BS
dx_bs=2*lambda;%BS inter-element distance
dx_user=2*lambda;%user inter-element distance
dx_ris=0.5*lambda;%ris inter-element distance
alpha=0.5;%reflection efficiency
Ropt = zeros(length(snrDbSet),numMC);
RpsInf= zeros(length(snrDbSet),numMC);
RpsInf1= zeros(length(snrDbSet),numMC);
%%
for mm = 1:numMC
    if mod(mm,10)==1
        mm
    end
    mm
    % generate chanel matrix
    Nc=3;Nray=5;
    G = channel_generation_ura(M,Mr,Nc,Nray,dx_ris,dx_user,lambda,p0,p1);
    H = channel_generation(M,Mt,d0,alpha,p0,lambda,dx_bs,dx_ris);
    Hd = zeros(Mr,Mt);%the direction channel between BS and User
    %% calculate rates
    for indxSnrDb = 1:length(snrDbSet)
        snrDb = snrDbSet(indxSnrDb);
        snrLin = db2pow(snrDb);
        %% upper bound
        Ropt(indxSnrDb,mm) = RIS_upperbound(G,H,snrLin,Ns,'rank1');
        %%  perfect PS based, No Switch H rank1
        RpsInf1(indxSnrDb,mm) = RIS_hrank1(G,H,snrLin);
        %% perfect PS based, No Switch
        RpsInf(indxSnrDb,mm) = RIS(G,H,snrLin,Ns,Hd,'psInf');
    end
end
%% det(Popt'*chanMat'*chanMat*Popt)
figure
width = 1.5;
plot(snrDbSet,mean(RpsInf,2),'b-v','LineWidth',width), hold on
plot(snrDbSet,mean(RpsInf1,2),'r-^','LineWidth',width), hold on
plot(snrDbSet,mean(Ropt,2),'k','LineWidth',width), hold on
legend('RIS scheme (b=\infty)','The power-like method in Section III.C','The Upper Bound','Location','NorthWest')
grid on
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bps/Hz)')
title(['URA in a Channel with ' num2str(Nc*Nray) ' Multipath Clusters, M_t=' num2str(Mt) ', M_r=' num2str(Mr) ', N_s=' num2str(Ns) ', M=' num2str(M) ', \beta=' num2str(alpha)]);
return

