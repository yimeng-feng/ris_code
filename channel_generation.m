
function [H1]=channel_generation(M,Mt,d0,alpha,m,lambda,dx_bs,dx_ris)
% ***************************************
% generate channel from bs to ris
%  author - Yimeng Feng
%  input- M: ris element number
%            Mt: bs element number
%            d0: distance between ris and bs
%            Nray: ray number
%            dx_t£ºinter-element distance for transmitter
%            dx_r£ºinter-element distance for receiver
%            lambda£ºwavelength
%            p0: line number for transmitter
%            p1: line number for receiver
%  output-H: Channel between BS to ris
%
%copyright - CSRL@Fudan,2021/01/18
%  ************************************
n=M/m;%the column of ris
H=zeros(M,Mt);
kk=1;
for ll=1:sqrt(Mt)
    for ii=1:sqrt(Mt)
        x1=-(sqrt(Mt)-1)/2*dx_bs+(ii-1)*dx_bs;
        y1=(sqrt(Mt)-1)/2*dx_bs-(ll-1)*dx_bs;
        jj=1;
        for mm=1:m
            for nn=1:n
                x2=-(n-1)/2*dx_ris+(nn-1)*dx_ris;
                y2=(m-1)/2*dx_ris-(mm-1)*dx_ris;
                d=sqrt((x2-x1)^2+(y2-y1)^2+d0^2);
                H(jj,kk)=exp(1i*2*pi*d/lambda)*sqrt(alpha/M);
                jj=jj+1;
            end
        end
        kk=kk+1;
    end
end
H1=H';
