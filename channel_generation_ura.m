function [H,alpha , Ar, At]=channel_generation_ura(Nt,Nr,Nc,Nray,dx_t,dx_r,lambda,p0,p1)
% ***************************************
%  generate URA channel
%  author - Yimeng Feng
%  input- Nt: transmitter antenna number
%            Nr: receiver antenna number
%            Nc: cluster number
%            Nray: ray number
%            dx_t: inter-element distance for transmitter
%            dx_r: inter-element distance for receiver
%            lambda: wavelength
%            p0: line number for transmitter
%            p1: line number for receiver
%  output-H: URA Channel 
%             alpha: multipath gain
%             Ar: antenna response for receiver
%             At: antenna response for transmitter
%copyright - CSRL@Fudan,2021/01/12
%  ************************************
angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
gamma = sqrt(1/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H
   
    for c = 1:Nc
        AoD_m = unifrnd(0,2*pi,1,2);
        AoA_m = unifrnd(0,2*pi,1,2);
        
        AoD(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(1),angle_sigma);
        AoD(2,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(2),angle_sigma);
        AoA(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoA_m(1),angle_sigma);
        AoA(2,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    end
    
    H = zeros(Nr,Nt);
    At = zeros(Nt,Nc*Nray);
    Ar = zeros(Nr,Nc*Nray);
    alpha=zeros(Nc*Nray,1);
    for j = 1:Nc*Nray
        At(:,j) = array_response(AoD(1,j),AoD(2,j),Nt,dx_t,lambda,p0); %UPA array response
        Ar(:,j) = array_response(AoA(1,j),AoA(2,j),Nr,dx_r,lambda,p1);
        alpha(j) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
        H = H + alpha(j) * Ar(:,j) * At(:,j)';
    end
    H = gamma * H;
    alpha=alpha*gamma;