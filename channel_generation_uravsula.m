function [H,Ar,alpha]=channel_generation_uravsula(Mt,Mr,Nc,Nray,dx_t,dx_r,lambda,p0)
% ***************************************
%  generate URA for ris and ula for bs and user' channel
%  author - Yimeng Feng
%  input- Mt: transmitter antenna number
%            Mr: receiver antenna number
%            Mr: receiver antenna number
%            dx_t: inter-element distance for transmitter
%            dx_r: inter-element distance for receiver
%            lambda: wavelength
%  output-H: Channel
%             Ar: antenna response for receiver
%             alpha: multipath gain
%copyright - CSRL@Fudan,2021/01/12
%  ************************************

angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
gamma = sqrt(1/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H

for c = 1:Nc
    AoD_m = unifrnd(0,2*pi,1,1);
    AoA_m = unifrnd(0,2*pi,1,2);
    
    AoD(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoA(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoA_m(1),angle_sigma);
    AoA(2,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoA_m(2),angle_sigma);
end

H = zeros(Mr,Mt);
At = zeros(Mt,Nc*Nray);
Ar = zeros(Mr,Nc*Nray);
alpha = zeros(Nc*Nray,1);
for j = 1:Nc*Nray
    At(:,j) = exp(1j*(2*pi*dx_t/lambda)*sin(AoD(1,j))*(0:Mt-1)'); %ULA array response
    Ar(:,j) = array_response(AoA(1,j),AoA(2,j),Mr,dx_r,lambda,p0); %UPA array response
    alpha(j) = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
    H = H + alpha(j) * Ar(:,j) * At(:,j)';
end
H = gamma * H;
end