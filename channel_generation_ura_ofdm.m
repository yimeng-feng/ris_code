function H=channel_generation_ura_ofdm(Nt,Nr,Nc,Nray,K,dx_t,dx_r,lambda,p0,p1)
% ***************************************
%  generate URA channel for ofdm system
%  author - Yimeng Feng
%  input- Nt: transmitter antenna number
%            Nr: receiver antenna number
%            Nc: cluster number
%            Nray: ray number
%            K: subcarrier number
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
Ht=zeros(Nr,Nt,Nc);
for c = 1:Nc
    AoD_m = unifrnd(0,2*pi,1,2);
    AoA_m = unifrnd(0,2*pi,1,2);
    
    AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
    AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
    AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    
    Ht(:,:,c) = zeros(Nr,Nt);
    At = zeros(Nt,Nc*Nray);
    Ar = zeros(Nr,Nc*Nray);
    for j = 1:Nray
        temp = (c-1)*Nray+j;
        At(:,temp) = array_response(AoD(1,j),AoD(2,j),Nt,dx_t,lambda,p0);
        Ar(:,temp) = array_response(AoA(1,j),AoA(2,j),Nr,dx_r,lambda,p1);
        alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
        Ht(:,:,c) = Ht(:,:,c) + alpha * Ar(:,temp) * At(:,temp)';
    end
end
H=zeros(Nr,Nt,K);
H1=zeros(Nt,Nr,K);
for k = 1:K
    H(:,:,k) = zeros(Nr,Nt);
    for c = 1:Nc
        H(:,:,k) = H(:,:,k) + Ht(:,:,c) * exp(-1i*2*pi/K*(k-1)*(c-1));
    end
    H(:,:,k) = H(:,:,k) * gamma;
    H1(:,:,k)=H(:,:,k)';
end