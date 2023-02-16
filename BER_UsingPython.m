close all
clear all
%% LoRa BER Underground Rayleigh Using Q function and Simplified Binomial
Pt=20;
BW=[125*10^3,250*10^3,500*10^3];
SF=[7,10,12];
NF=6;
mv=[0.05,0.1,0.15];
S=0.9;
C=0.1;
%% Load data processed by python
BER_Diversity_Using_Python_SF7=importdata('./BER_Rayleigh_sf7.txt');
BER_Diversity_Using_Python_SF7=erase(BER_Diversity_Using_Python_SF7,["[","]"]);
BER_Diversity_Using_Python_SF7=split(BER_Diversity_Using_Python_SF7,',');
BER_Diversity_Using_Python_SF7=str2num(char(BER_Diversity_Using_Python_SF7));
BER_Diversity_Using_Python_SF10=importdata('./BER_Rayleigh_sf10.txt');
BER_Diversity_Using_Python_SF10=erase(BER_Diversity_Using_Python_SF10,["[","]"]);
BER_Diversity_Using_Python_SF10=split(BER_Diversity_Using_Python_SF10,',');
BER_Diversity_Using_Python_SF10=str2num(char(BER_Diversity_Using_Python_SF10));
BER_Diversity_Using_Python_SF12=importdata('./BER_Rayleigh_sf12.txt');
BER_Diversity_Using_Python_SF12=erase(BER_Diversity_Using_Python_SF12,["[","]"]);
BER_Diversity_Using_Python_SF12=split(BER_Diversity_Using_Python_SF12,',');
BER_Diversity_Using_Python_SF12=str2num(char(BER_Diversity_Using_Python_SF12));

for k=1:1:140
Distance(k)=k*0.1;
SNR=UndergroundPassLoss(Pt,Distance(k),mv(3),S,C,BW(1),NF);
% Pb_Rayleigh_BinomialSimplified_SF7(k)=Pb_Rayleigh_Using_Simplified_Binomial(SF(1),SNR);
% Pb_Rayleigh_BinomialSimplified_SF10(k)=Pb_Rayleigh_Using_Simplified_Binomial(SF(2),SNR);
% Pb_Rayleigh_BinomialSimplified_SF12(k)=Pb_Rayleigh_Using_Simplified_Binomial(SF(3),SNR);
Pb_Rayleigh_Q_Underground_SF7(k)=Pb_Rayleigh_Using_Q_Fuction(SF(1),SNR);
Pb_Rayleigh_Q_Underground_SF10(k)=Pb_Rayleigh_Using_Q_Fuction(SF(2),SNR);
Pb_Rayleigh_Q_Underground_SF12(k)=Pb_Rayleigh_Using_Q_Fuction(SF(3),SNR);
% Pb_Rayleigh_Asymptotic_SF12(k)=Pb_Rayleigh_AsymptoticBER(SF(3),SNR);
% Pb_Rayleigh_Binomial(k)=Pb_Rayleigh_Using_Binomial(SF(1),SNR);
% Pb_Rayleigh_Diversity_Binomial(k)=Pb_Rayleigh_Diversity_Using_Binomial(SF(1),SNR,SNR,SNR);
end


figure
semilogy(Distance,Pb_Rayleigh_Q_Underground_SF7,'-r')
hold on
semilogy(Distance,Pb_Rayleigh_Q_Underground_SF10,'-c')
hold on
semilogy(Distance,Pb_Rayleigh_Q_Underground_SF12,'-b')
% hold on
% semilogy(Distance,Pb_Rayleigh_Binomial,'--g')
% hold on
% semilogy(Distance,Pb_Rayleigh_Asymptotic_SF12,'-k')
hold on
semilogy(Distance,BER_Diversity_Using_Python_SF7,'-.r')
hold on
semilogy(Distance,BER_Diversity_Using_Python_SF10,'-.c')
hold on
semilogy(Distance,BER_Diversity_Using_Python_SF12,'-.b')

% legend('SF=7','SF=8','SF=9','SF=10','SF=11','SF=12','BinomialSimplified')



function  SNR=UndergroundPassLoss(Pt,d,mv,S,C,BW,NF)
    Gt=0;
    Gr=0;
    f=900e6;
    ew_inf=4.9;
    ew0=80.1;
    e_air=8.854e-12;
    rho_s=2.65;% per centimeter
    rho_b=1.5;
    delta_eff=0.0467+0.2204*rho_b-0.4111*S+0.6614*C;
    e_water_real=ew_inf+(ew0-ew_inf)/(1+(0.58*10^(-10)*f)^2);
    e_water_imag=0.58*10^(-10)*f*(ew0-ew_inf)/(1+(0.58*10^(-10)*f)^2)+delta_eff/(2*pi*f*e_air)*(rho_s-rho_b)/(rho_s*mv);
    e0=1.15*(1+rho_b/rho_s*(((1.01+0.44*2.65)^2-0.062)^0.65-1)+mv^(1.2748-0.519*S-0.152*C)*e_water_real^0.65-mv)^(1/0.65)-0.68;
    e1=((mv^(1.33797-0.603*S-0.166*C))*(e_water_imag)^0.65)^(1/0.65);
    e0=e0*e_air;
    e1=e1*e_air;
    mu=1.0006*4*pi*1e-7;
    alpha=2*pi*f*sqrt((mu*e0/2)*(sqrt(1+(e1/e0)^2)-1));
    beta=2*pi*f*sqrt((mu*e0/2)*(sqrt(1+(e1/e0)^2)+1));
    lamda=2*pi/beta;
    Lug=6.4+20*log10(d)+20*log10(beta)+8.69*alpha*(d);
    Pr=Pt+Gt+Gr-Lug;
    SNR=Pr+174-10*log10(BW)-NF;
end


function Pb_Rayleigh_Binomial=Pb_Rayleigh_Diversity_Using_Binomial(SF,Pr1,Pr2,Pr3)
% SNR_Underground_dB1=Pr1+174-10*log10(BW)-NF;
SNR_Underground1=10^(Pr1/10);
SNR1=SNR_Underground1*2^SF;
% SNR_Underground_dB2=Pr2+174-10*log10(BW)-NF;
SNR_Underground2=10^(Pr2/10);
SNR2=SNR_Underground2*2^SF;
% SNR_Underground_dB3=Pr3+174-10*log10(BW)-NF;
SNR_Underground3=10^(Pr3/10);
SNR3=SNR_Underground3*2^SF;
temp=0;
for q=1:1:2^SF-1
Cnk=nchoosek(2^SF-1,q);
c2=(1+q/(q+1)*SNR1)^(-1)*(1+q/(q+1)*SNR2)^(-1)*(1+q/(q+1)*SNR3)^(-1);
temp=temp+(-1)^(q+1)*Cnk/(q+1)*c2;
end
Pb_Rayleigh_Binomial=(2^(SF-1)/(2^(SF)-1))*temp;
end

function Pb_Rayleigh_Binomial=Pb_Rayleigh_Using_Binomial(SF,SNR)
SNR_Underground1=10^(SNR/10);
SNR1=SNR_Underground1*2^SF;
temp=0;
for q=1:1:2^SF-1
Cnk=nchoosek(2^SF-1,q);
c2=(1+q/(q+1)*SNR1)^(-1);
temp=temp+(-1)^(q+1)*Cnk/(q+1)*c2;
end
Pb_Rayleigh_Binomial=(2^(SF-1)/(2^(SF)-1))*temp;
end

function Pb_Rayleigh_Asymptotic=Pb_Rayleigh_AsymptoticBER(SF,Pr1)
% SNR_Underground_dB1=Pr1+174-10*log10(BW)-NF;
SNR_Underground1=10^(Pr1/10);
SNR1=SNR_Underground1*2^SF;
% SNR_Underground_dB2=Pr2+174-10*log10(BW)-NF;
% SNR_Underground2=10^(Pr2/10);
% SNR2=SNR_Underground2*2^SF;
% SNR_Underground_dB3=Pr3+174-10*log10(BW)-NF;
% SNR_Underground3=10^(Pr3/10);
% SNR3=SNR_Underground3*2^SF;
gamma=0.577;
Pb_Rayleigh_Asymptotic=(2^(SF-1)/(2^(SF)-1))*(gamma+log(2^(SF)-1))/(SNR1+1);%/(SNR2*2^SF+1)/(SNR3*2^SF+1);
end
function Pb_Rayleigh_Binomial=Pb_Rayleigh_Using_Simplified_Binomial(SF,SNR)
SNR_Underground=10^(SNR/10);
SNR1=SNR_Underground;
SNR1_gamma1=gammaln(2^(SF));
SNR1_gamma2=gammaln((2+SNR1*2^SF)/(1+SNR1*2^SF));
SNR1_gamma3=gammaln((1+2^(SF)+2^(SF)*SNR1*2^SF)/(1+SNR1*2^SF));
temp1=exp(SNR1_gamma1+SNR1_gamma2-SNR1_gamma3);
Pb_Rayleigh_Binomial=(2^(SF-1)/(2^(SF)-1))*(1-temp1);
end

function Pb_Rayleigh_Q_Underground=Pb_Rayleigh_Using_Q_Fuction(SF,SNR)
%     SNR_Underground_dB=Pr+174-10*log10(BW)-NF;
    SNR_Underground=10^(SNR/10);
    gamma_eff_Underground=SNR_Underground*2^SF;
    H=0;
    for j=1:1:2^SF-1
    H=H+1/j;
    end
    part1_Underground=qfunc(-sqrt(2*H));
    part2_Underground=sqrt(gamma_eff_Underground/(gamma_eff_Underground+1)*exp(-2*H/(2*(gamma_eff_Underground+1))));
    part3_Underground=qfunc(sqrt((gamma_eff_Underground+1)/(gamma_eff_Underground))*(-sqrt(2*H)+(sqrt(2*H)/(gamma_eff_Underground+1))));
    Pb_Rayleigh_Q_Underground=0.5*(part1_Underground-part2_Underground*part3_Underground);
end