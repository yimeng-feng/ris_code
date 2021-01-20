function [C,P,ldpow] = calc_capacity(H,Ns,snrLin,useWf)
% ***************************************
% calculate the mimo capacity log|I+HPP'H'|
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
if nargin<4
    useWf = 1;
end
[~,S,V] = svd(H);
sv = diag(S);
sv = sv(1:Ns);
if useWf
    ldpow = calc_waterfilling(sv,snrLin,1);
else
    ldpow = snrLin/Ns*ones(Ns,1);
end
P = V(:,1:Ns)*diag(sqrt(ldpow));
HP = H*P;
C = real(log2(det(eye(Ns)+HP'*HP)));