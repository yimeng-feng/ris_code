function ldpow=calc_waterfilling_ofdm(sv,snrLin,Ns,N,sigma2)
% ***************************************
% calculate waterfilling power
%  author - Yimeng Feng
%  input- sv: singular values of the channel
%            snrLin: input power one subcarrier
%            Ns: data stream
%            N: subcarrier number
%            sigma2: noise power
%  output- ldpow: output power
%
%copyright - CSRL@Fudan,2021/01/18

sv2=sv.^2;
b=sort(sv2(:),'descend');
p=snrLin*N;
for nn=1:N
    ind(:,nn)=find(ismember(b,sv2(:,nn)));
end
for ii=N*Ns:-1:1
    mu = p/ii+sum(sigma2./b(1:ii))/ii;
    if mu> sigma2/b(ii)
        break
    end
end
ldpow = zeros(Ns,N);
b(ii+1:end)=0;
for n = 1:N
    ind_n = ind(:,n);
    len = length(ind_n(ind_n<ii+1));
    ldpow(1:len,n) = mu-sigma2./b(ind(1:len,n));  
end

end