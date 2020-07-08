clear all
close all
clc

h=6.626e-34; %[J*s]
c=3e8; %[m/s]
qe=1.6e-19; %[C]
Eg=1.95*qe; %[J]
lm=550; %[nm]
FWHM=300; %[nm]
sigma=FWHM/sqrt(2*log(2)); %[nm]
lg=1:1:1500; %[nm]
s2=sigma*sqrt(2);
u=(lm*(erf(lm/s2)+erf((lg-lm)/s2))+s2/sqrt(pi)*(exp(-(lm/s2)^2)-exp(-((lg-lm)/s2).^2)))./(lg*(1+erf(lm/s2)));

figure
plot(lg,u)
xlabel('\lambda [nm]')
ylabel('u')
title('PCE (at the ideal Shockley-Queisser limit) vs. energy gap wavelength')
lideal = find(max(u)==u); %Wavelength for ideal Shockley-Queisser limit [nm]
hold all
line([lideal lideal],[0 u(lideal)],'Color','r')
%Desired bandgap for the maximum PCE would be Egideal
Egideal = h*c/lideal/qe*1e9 %[eV]
