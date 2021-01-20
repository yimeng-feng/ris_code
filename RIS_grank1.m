function Ropt= RIS_grank1(G,H,snrLin)
% ***************************************
% algorithm for maximizing RIS phases when G rank1
%  author - Yimeng Feng
%  input-  G: channel between ris and user
%             H: channel between bs and ris
%             snrLin: signal noise ratio
%  output- Ropt: mimo capacity
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
maxIterNum = 1000;
f = zeros(maxIterNum,1);
[U,~,V] = svd(G);
Atheta = U(:,1)./abs(U(:,1));
Aphi = V(:,1)./abs(V(:,1));
A = conj(H)*diag(Aphi');
B = A'*A;
[U1,~,~] = svds(B,1);
d = exp(1j*angle(U1(:,1)));%initialization for D
f(1) = real(d'*B*d);
ii=1;
pre = 0;
while abs(f(ii)-pre)>1e-4&&ii<maxIterNum
    pre=f(ii);
    d=exp(1j*angle(B*d));
    ii=ii+1;
    f(ii)=real(d'*B*d);
end
p1 = H*diag(d)'*Aphi;
p = p1/norm(p1)*sqrt(snrLin);%popt
% figure
% plot(f)
rho=norm(G*diag(d)*H'*p)^2;

Ropt=log2(1+rho);
end