function Ropt=RIS_upperbound(G,H,snrLin,Ns,type)
% ***************************************
% calculate the upperbound
%  author - Yimeng Feng
%  input- G: Channel between ris and user
%            H: Channel between bs and ris
%            snrLin: signal noise ratio
%            Ns: data stream
%            type: rank1 case for channel H or G rank 1
%  output-Ropt: capacity upperbound
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
if nargin<5
    type = 'full rank';
end
switch type
    case 'rank1'
        [~,S1,~]=svd(G*G');
        s1=max(diag(S1));
        [~,S2,~]=svd(H*H');
        s2=max(diag(S2));
        Ropt=log2(1+s1*s2*snrLin);
    case'full rank'
        [~,~,V1]=svd(G);
        [U2,~,~]=svd(H');
        U=V1*U2';
        %calculate P by waterfilling
        Ropt=calc_capacity(G*U*H',Ns,snrLin,1);
end

end
