function Ropt=RIS_ofdm(G,H,snrLin,Ns,Hd,type)
% ***************************************
% algorithm for maximizing RIS phases in ofdm system
%  author - Yimeng Feng
%  input-  G: channel between ris and user
%             H: channel between bs and ris
%             snrLin: signal noise ratio
%             Ns: data stream
%             Hd: channel between bs and user
%             type: psInf infinite resolution%
%                      ps2 2 bit resolution
%                      ps1 1bit resolution
%  output- Ropt: mimo capacity
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
%tic
numlter=1000;
R=zeros(numlter,1);
[Mr,M,N] = size(G);
[Mt,~,~]=size(H);
d = randn(M,1)+1i*randn(M,1);%initialization
D = diag(d./abs(d));
H1 = zeros(Ns,M,N);
C = zeros(Mr,Ns,N);
A = zeros(Mr,Mr,N);
invA = zeros(Mr,Mr,N);
HH1 = zeros(Mr,Mt,N);
for nn=1:N
    HH1(:,:,nn)=G(:,:,nn)*D*H(:,:,nn)'+Hd(:,:,nn);
end
[R(1),P]=calc_capacity_ofdm(HH1,Ns,snrLin,1);%Calculate p and capacity

for nn=1:N
    H1(:,:,nn)=P(:,:,nn)'*H(:,:,nn);
end
R_pre=0; ii=1;
while abs(R(ii)-R_pre)>=1e-3 && ii<=20
    %fix FBB, calculate D
    R1=zeros(M,1);
    ii=ii+1;
    for nn=1:N
        C(:,:,nn)=G(:,:,nn)*D*H1(:,:,nn)'+Hd(:,:,nn)*P(:,:,nn);
        A(:,:,nn)=eye(Mr)+C(:,:,nn)*C(:,:,nn)';
        [U,S,~]=svd(A(:,:,nn));
        sv=diag(S);
        invA(:,:,nn)=U*diag(1./sv)*U';%the inversion of A
    end
    for mm=1:M
        R2=zeros(1,100);
        e1=zeros(Mr,N);
        e2=zeros(Mr,N);
        e3=zeros(N,1);
        e4=zeros(Mr,N);
        e5=zeros(N,1);
        e6=zeros(N,1);
        c1=zeros(N,1);
        c2=zeros(N,1);
        gn=zeros(N,1);
        cn=zeros(N,1);
        for nn=1:N
            G1=G(:,:,nn);
            HH=H1(:,:,nn);
            e1(:,nn)=invA(:,:,nn)*G1(:,mm);
            e2(:,nn)=C(:,:,nn)*HH(:,mm);
            e3(nn)=HH(:,mm)'*HH(:,mm);
            e4(:,nn)=invA(:,:,nn)*e2(:,nn);
            e5(nn)=e2(:,nn)'*e4(:,nn);
            e6(nn)=G1(:,mm)'*e1(:,nn);
            c1(nn)=e1(:,nn)'*e2(:,nn);
            c2(nn)=c1(nn)*c1(nn)'+e6(nn)*(e3(nn)-e5(nn));
            gn(nn)=c1(nn)-c2(nn)*D(mm,mm);
            cn(nn)=2*c2(nn)+1-2*real(D(mm,mm)'*c1(nn));
        end
        theta=angle(D(mm,mm));
        %% optimize theta by back-tracking line search
        R2(1)=15;R2_pre=0.5;
        jj=1;
        while abs(R2(jj)-R2_pre)>=1e-3&& jj<=10
            jj=jj+1;
            df1=0;
            f0=0;
            ff=0;
            for nn=1:N
                dfa=-1j*exp(-1j*theta)*gn(nn)+1j*exp(1j*theta)*gn(nn)';
                dfb=real(cn(nn)+2*real(gn(nn)*exp(-1j*theta)));
                df1=df1-dfa/dfb;%first-order derivative
                f0=f0-log(dfb);
                ff=ff-log2(dfb);
            end
            R2(jj)=ff;
            R2_pre=R2(jj-1);
            %back tracking line search
            s=1;b1=0.5;a1=0.3;
            v=-df1;
            theta1=theta+s*v;
            f1=0;
            for nn=1:N
                dfb=real(cn(nn)+2*real(gn(nn)*exp(-1j*theta1)));
                f1=f1-log(dfb);
            end
            while f1>f0+a1*s*df1*v
                s=b1*s;
                theta1=theta+s*v;
                f1=0;
                for nn=1:N
                    f1=f1-log(cn(nn)+2*real(gn(nn)*exp(-1j*theta1)));
                end
                if s<=1e-7
                    break
                end
            end
            theta=theta+s*v;
        end
        %% update theta
        t1=exp(1j*theta)-D(mm,mm);
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
        e7 = t1*t1';
        e8=zeros(N,1);
        for nn=1:N
            G1=G(:,:,nn);
            HH=H1(:,:,nn);
            e8(nn)=1/2*e7*e3(nn);
            % matrix inversion lemma
            Af= t1'*e4(:,nn)+ e8(nn)*e1(:,nn);
            gAf =e8(nn)*e6(nn)+t1'*c1(nn);
            fAf=e7*e5(nn)+t1*e8(nn)*c1(nn)'+e8(nn)*gAf;
            a=real(e6(nn));
            b=real(fAf);
            c=gAf+1;
            d=1/(a*b-c*c');
            t11=d*b;
            t12=-d*c;
            t22=d*a;
            TAf=t12*e1(:,nn)*Af';
            Ti=t11*(e1(:,nn)*e1(:,nn)')+t22*(Af*Af')+TAf+TAf';
            invA(:,:,nn)=invA(:,:,nn)-Ti;
            C(:,:,nn)=C(:,:,nn)+t1*G1(:,mm)*HH(:,mm)';%update C
            HP=(G(:,:,nn)*D*H(:,:,nn)'+Hd(:,:,nn))*P(:,:,nn);
            R1(mm)=R1(mm)+1/N*real(log2(det(eye(Mr)+HP*HP')));
        end
    end
    %     figure,plot(R1(1:mm-1))
    HH1=zeros(Mr,Mt,N);
    %% fix D, calculate p by waterfilling
    for nn=1:N
        HH1(:,:,nn)=G(:,:,nn)*D*H(:,:,nn)'+Hd(:,:,nn);
    end
    R_pre=R(ii-1);
    [R(ii),P]=calc_capacity_ofdm(HH1,Ns,snrLin,1);
    for nn=1:N
        H1(:,:,nn)=P(:,:,nn)'*H(:,:,nn);
    end
end
% figure,plot(R(1:ii))
Ropt=R(ii);
%toc
end


