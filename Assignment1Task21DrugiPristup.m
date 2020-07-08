clear all
close all
clc

h=6.626e-34; %[J*s]
c=299792458; %[m/s]
qe=1.6e-19; %[C]
kB=1.38e-23; %[J/K]
Eg=1.95*qe; %[J]
lm=550; %[nm]
FWHM=300; %[nm]
A=1.5; %[W/(m^2 nm)]
sigma=FWHM/sqrt(8*log(2)); %[nm]
lg=1:1:3000; %[nm]
s2=sigma*sqrt(2);
u=(lm*(erf(lm/s2)+erf((lg-lm)/s2))+s2/sqrt(pi)*(exp(-(lm/s2)^2)-exp(-((lg-lm)/s2).^2)))./(lg*(1+erf(lm/s2)));

%% Plotting the spectral irradiance reaching the planet surface
lambda = 1:1:2*lm;
elambda = A*exp(-(lambda-lm).^2/(2*sigma^2)); %Gaussian spectrum of NUS EM radiation reaching planet OTATOP surface
T = 5800; %[K], black body effective temperature of SUN 
el = 2*pi*h*c^2./((lg*(1e-9)).^5.*(exp(h*c./(lg*(1e-9)*kB*T))-1)); %Emission spectrum of SUN EM radiation
figure(1)
R_E = 6371; %[km]
R_ES = 150e6; %[km]
R_S = 696000; %[km]
koef = sum(el)*(1e-9)/(5.67*1e-8)/5778^4
plot(lg,el,'LineWidth',1.5);
xlabel('\lambda [nm]')
ylabel('e_{\lambda}(\lambda) [W/(m^2 nm)]')
title('Spectral emissive power of the Sun')

figure(2)
plot(lambda,elambda,'LineWidth',1.5)
xlim([0 2*lm])
xlabel('\lambda [nm]')
ylabel('I_{\lambda}(\lambda) [W/(m^2 nm)]')
title('Spectral irradiance reaching the planet surface')

%% Plotting ultimate efficiency vs lambda_gap
figure(3)
plot(lg,u,'LineWidth',1.5)
xlabel('\lambda [nm]')
ylabel('u')
title('PCE (at the ideal Shockley-Queisser limit) vs energy gap wavelength')
lideal = find(max(u)==u); %Wavelength for ideal Shockley-Queisser limit [nm]
hold all
line([lideal lideal],[0 u(lideal)],'Color','r')
line([0 lideal],[u(lideal) u(lideal)],'Color','r')
%Desired bandgap for the maximum PCE would be Egideal
Egideal = h*c/lideal/qe*1e9 %[eV]
PCEideal = max(u)
PCE637 = u(637)
%% Power conversion efficiency at lambda_gap=h*c/E_gap
FF = 0.7;
EQE = 0.85;
lg = h*c/Eg*1e9;
u=(lm*(erf(lm/s2)+erf((lg-lm)/s2))+s2/sqrt(pi)*(exp(-(lm/s2)^2)-exp(-((lg-lm)/s2).^2)))./(lg*(1+erf(lm/s2)));
etaabs=(lm*(erf(lm/s2)+erf((lg-lm)/s2))+s2/sqrt(pi)*(exp(-(lm/s2)^2)-exp(-((lg-lm)/s2).^2)))./(lm*(1+erf(lm/s2))+s2/sqrt(pi)*exp(-(lm/s2)^2));
IQE = EQE/etaabs;
PCE1 = FF*EQE*u;
PCE = FF*EQE*(lm*(erf(lm/s2)+1)+s2/sqrt(pi)*(exp(-(lm/s2)^2)))./(lg*(1+erf(lm/s2)))