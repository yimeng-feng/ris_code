function [C,P,ldpow] = calc_capacity_ofdm(H,Ns,snrLin,useWf)
% ***************************************
% calculate the sum mimo capacity for ofdm system
%  author - Yimeng Feng
%  input- H: channel H
%            Ns: data stream
%            snrLin: signal noise ratio
%            useWf: 1 waterfilling
%                        0 no waterfilling
%  output-C: mimo capacity
%              P: digital precoder
%              Idpow: power for each stream
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
C=0;
[~,Mt,N]=size(H);
P=zeros(Mt,Ns,N);
V=zeros(Mt,Mt,N);
for nn=1:N
    [~,S,V(:,:,nn)] = svd(H(:,:,nn));
    s=diag(S);
    sv(:,nn)=s(1:Ns);
end
if useWf
    ldpow = calc_waterfilling_ofdm(sv,snrLin,Ns,N,1);
else
    ldpow = snrLin/Ns*ones(Ns,N);
end
for nn=1:N
    P(:,:,nn) = V(:,1:Ns,nn)*diag(sqrt(ldpow(:,nn)));
    HP = H(:,:,nn)*P(:,:,nn);
    C =C+ 1/N*real(log2(det(eye(Ns)+HP'*HP)));%calculate rate
end
end