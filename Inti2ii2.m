%% vicinity collision
% version 2.0 
% this version provides an analytical formula
function I_T=Inti2ii2(T,Z_M,LD)
% inputs£ºquantity, unit
% T: temperature, K
% Z_M: average ionization, 
% LD, Debye length, m

% -----output-----£º
% I_T£º integral term, m^2

%% physical constant
epsilon0=8.8542e-12; % F/m, electric constant in vacuum
elect=1.6022e-19; % C, elementary charge
kb=1.3807e-23; % J/K, Boltzmann constant
%% q^2, charge square
% Z_M=Z_M*0.02;
qs=(Z_M*elect).^2/4/pi/epsilon0;
% q=sqrt(qs);

%% parameters of the hyperbola
% equation: x^2/a^2-y^2/b^2=1
r0=LD; % m, cut off distance
a=qs./(3*kb.*T+2.*qs./r0); % m, semi-major axis

%% ingegral results
K=((r0-a)./a).^2;
I_T=2.*r0.^2./(1+sqrt(K)).^2.*K./(K-1).^2.*(  2.*(1-K)+(1+K).*log(K)    );


end




