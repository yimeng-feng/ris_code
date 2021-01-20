%%
% currentFolder = pwd;
% addpath(genpath(currentFolder));
%%
clc, clear
Mt = 64;%BS antenna
Mr =64;%user antenna
Ns=4;%data stream
M =256;%ris element number
numMC =10;
snrDbSet =-20:5:10;
p0=16;
c=3e8;%speed of light
fc=7e10;%carrier frequency
lambda=c/fc;%wavelength
d0=1;%distance between ris and BS
dx_bs=2*lambda;%BS inter-element distance
dx_user=2*lambda;%user inter-element distance
dx_ris=0.5*lambda;%ris inter-element distance
alpha=0.5;%reflection efficiency
Ropt = zeros(length(snrDbSet),numMC);
RpsInf= zeros(length(snrDbSet),numMC);
Rmani= zeros(length(snrDbSet),numMC);
%%
for mm = 1:numMC
    if mod(mm,10)==1
        mm
    end
    mm
    % generate chanel matrix
    Nc=3;%cluster number
    Nray=5;%ray number
    [H0,Ar,ar]=channel_generation_uravsula(Mt,M,Nc,Nray,dx_bs,dx_ris,lambda,p0);
    H=H0';
    [G0,At,at]=channel_generation_uravsula(Mr,M,Nc,Nray,dx_user,dx_ris,lambda,p0);
    G=G0';
    Hd=zeros(Mr,Mt);
    %% calculate rates
    for indxSnrDb = 1:length(snrDbSet)
        snrDb = snrDbSet(indxSnrDb);
        snrLin = db2pow(snrDb);
        %% upper bound
        Ropt(indxSnrDb,mm) = RIS_upperbound(G,H,snrLin,Ns);
        %% perfect PS based, No Switch
        RpsInf(indxSnrDb,mm) = RIS(G,H,snrLin,Ns,Hd,'psInf',1);
        %% manicode for t-svd
        trans_Pt=snrLin;
        sigma2=1;
        [H_man,v_man] =  test_Mani_0714( G,H',M, Ns,Ar,ar,at,At,trans_Pt,sigma2);
        D=diag(v_man');
        Rmani(indxSnrDb,mm)=calc_capacity(G*D*H',Ns,snrLin,1);
    end
end
%% det(Popt'*chanMat'*chanMat*Popt)
figure
width = 1.5;
plot(snrDbSet,mean(Rmani,2),'r-s','LineWidth',width), hold on
plot(snrDbSet,mean(RpsInf,2),'b-d','LineWidth',width), hold on
plot(snrDbSet,mean(Ropt,2),'k','LineWidth',width), hold on
legend('T-SVD in [14]','Our proposed algorithm (b=\infty)','The Upper Bound','Location','NorthWest')
grid on
xlabel('SNR (dB)')
ylabel('Spectral Efficiency (bps/Hz)')
title(['URA in a Channel with ' num2str(Nc*Nray) ' Multipath Clusters, M_t=' num2str(Mt) ', M_r=' num2str(Mr) ', N_s=' num2str(Ns) ', M=' num2str(M) ', \beta=' num2str(alpha)]);
return
