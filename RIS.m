function Ropt=RIS(G,H,snrLin,Ns,Hd,type,up1rank)
% ***************************************
% algorithm for maximizing RIS phases
%  author - Yimeng Feng
%  input-  G: channel between ris and user
%             H: channel between bs and ris
%             snrLin: signal noise ratio
%             Ns: data stream
%             Hd: channel between bs and user
%             type: psInf infinite resolution%
%                      ps2 2 bit resolution
%                      ps1 1bit resolution
%             up1rank: 1 use 1 rank for update invA
%  output- Ropt: mimo capacity
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
if nargin<7
    up1rank=0;
end
numlter=1000;
R=zeros(numlter,1);
[Mr,M] = size(G);
%% initialization
d = randn(M,1)+1i*randn(M,1);
D = diag(d./abs(d));
HP = G*D*H'+Hd;
[ R(1), P ] = calc_capacity(HP,Ns,snrLin,1);% initial point
R_pre = 0; ii = 1;
while abs(R(ii)-R_pre)>=1e-4 && ii<=15
    ii=ii+1;
    %fix FBB, calculate D
    H1=P'*H;
    %         R1=zeros(M,1);
    C=zeros(Mr,Ns);
    for jj=1:M
        C=C+G(:,jj)*D(jj,jj)*H1(:,jj)';
    end
    A=eye(Mr)+C*C';
    invA=inv(A);%the inversion of A
    for mm = 1:M
        %e1~e8 are Intermediate variables to decrease computation complexity
        e1 = invA*G(:,mm);
        e2 = C*H1(:,mm);
        e3 = H1(:,mm)'*H1(:,mm);
        e4 = invA*e2;
        e5 = real(e2'*e4);
        e6 = real(G(:,mm)'*e1);
        c1 = e1'*e2;
        c2 = real(c1*c1'+e6*(e3-e5));
        theta = angle(c1-D(mm,mm)*c2);
        t1 = exp(1j*theta)-D(mm,mm);
        %update D(mm,mm)
        switch type
            case 'psInf'
                D(mm,mm)=exp(1j*theta);
            case 'ps2'
                D(mm,mm)=exp(1j*(2*pi*round((theta*2^2)/(2*pi)))/2^2);
            case'ps1'
                D(mm,mm)=exp(1j*(2*pi*round((theta*2^1)/(2*pi)))/2^1);
        end
        if mm==M
            break;
        end
        %% update invertion of A(invA) and C
        e7 = real(t1*t1');
        e8 = real(1/2*e7*e3);
        if up1rank==1
            % matrix inversion lemma by rank1 update
            p=t1'*C*H1(:,mm)+1/2*(t1*t1')*G(:,mm)*(H1(:,mm)'*H1(:,mm));
            invD=invA-(invA*G(:,mm)*p'*invA)/(1+p'*invA*G(:,mm));
            invA=invD-(invD*p*G(:,mm)'*invD)/(1+G(:,mm)'*invD*p);
        else
            % matrix inversion lemma
            Af= t1'*e4+ e8*e1;
            gAf =e8*e6+t1'*c1;
            fAf=e7*e5+t1*e8*c1'+e8*gAf;
            a=real(e6);
            b=real(fAf);
            c=gAf+1;
            d=1/(a*b-c*c');
            t11=d*b;
            t12=-d*c;
            t22=d*a;
            TAf=t12*e1*Af';
            Ti=t11*(e1*e1')+t22*(Af*Af')+TAf+TAf';
            invA=invA-Ti;
        end
        
        C=C+t1*G(:,mm)*H1(:,mm)';%update C
        %         HP=G*D*H'*P;
        %         R1(mm)=real(log2(det(eye(Mr)+HP*HP')));
    end
    %        figure,plot(R1(1:mm-1))
    %% fix D, calculate p by waterfilling
    HP=G*D*H'+Hd;
    [R(ii),P]=calc_capacity(HP,Ns,snrLin,1);
    R_pre=R(ii-1);
end
%       figure(5),plot(R(1:ii),'b-d','LineWidth',1.5),hold on
%       xlabel('iteration number')
%       ylabel('Spectral Efficiency (bps/Hz)')
%       title(['URA in a Channel with 15 Multipath Clusters, M_t=64, M_r=' num2str(Mr) ', N_s=' num2str(Ns)])
Ropt=R(ii);
end

