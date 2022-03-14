%% Calculate the average ionization employing TF model
% the fitting parameters are from 
% Chinese Journal of Computational Physics Vol.34, No.5, PP505, 2017
% on Page 511
% inputs: quantity, unit
% Z: nuclear charge, 
% Am: atomic weight, 
% n: number density of atom, m^-3
% T: temperature, K

% ----- output -----:
% Z_M: average ionization
%% function
function Z_M = Z_TF(Z,Am,n,T)

Na=6.02214e23; % /mol, Avogadro number 
elect=1.6022e-19; % C, elementary charge
kb=1.3807e-23; % J/K, Boltzmann constant
Te=T*kb/elect; % eV, temperature
Z_S=Z^(4/3); % Z_star
rho=n*Am*1e-3/Na*1e-3; %g/cm^3, mass density
%% fitting parameters

alpha=14.3139;
beta=0.6624;

a1=3.323e-3;
a2=0.97183;
a3=9.261e-5;
a4=3.101;

b1=-1.763;
b2=1.4317;
b3=0.3154;

c1=-0.366667;
c2=0.9833;

%% calculated parameters
A=a1.*(Te/Z_S).^a2+a3.*(Te/Z_S).^a4;
B=-exp(b1+b2.*(Te./(Te+Z_S))+b3.*(Te./(Te+Z_S)).^7);
C=c1.*(Te./(Te+Z_S))+c2;
x=alpha.*((rho/Z/Am).^C+A.^C.*(rho/Z/Am).^(B.*C)).^(beta./C);

%% Z_mean - T
Z_M=x*Z./(1+x+sqrt(1+x*2));

end

