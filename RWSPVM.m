% --------------Copyright(c)--------------
% File name: RWSPVM.m
% Version: 3.3
% Created by: Yuqing Cheng, ORCID: 0000-0002-3371-9164
% Created date: March 10, 2020
% Last modified date: January 10, 2022
% Descriptions: This is the program to calculate the viscosity in warm
% temperature range, employing our random-walk shielding-potential
% viscosity model (RWSP-VM). We give 4 examples: Al, Fe, U and Be. Other
% metals are also valid as long as one inputs their basic quantities.
% Moreover, the temperature range can be changed according to the requirements. 
% -----------------------------------------

%% Random-Walk Shielding-Potential Vicosity Model

%% RWSP-VM
% calculating viscosity of the specified metal in a flash
clc;clear; % clear the screen and the parameters.

%% Physical constants
epsilon0=8.8542e-12; % F/m, electric constant in vacuum
elect=1.6022e-19; % C, elementary charge
kb=1.3807e-23; % J/K, Boltzmann constant
Na=6.02214e23; % mol^-1, Avogadro number 
me=9.10938215e-31; % kg, mass of electron
h=6.62606896e-34; % J.s, Planck constant
%% Al
Z=13; % nuclear charge number
Am=26.982; % atomic weight
M=Am*1e-3; % kg/mol, molar mass
m=M/Na; % kg, mass of Al atom

Rho=2700; % kg/m^3, mass density
%% Fe
% Z=26; % nuclear charge number
% Am=55.845; % atomic weight
% M=Am*1e-3; % kg/mol, molar mass
% m=M/Na; % kg, mass of the atom
% 
% Rho=7.9e3;% kg/m^3, mass density

%% U
% Z=92; % nuclear charge number
% Am=238; % atomic weight
% M=Am*1e-3; % kg/mol, molar mass
% m=M/Na; % kg, mass of the atom
% 
% Rho=94.65e3; % kg/m^3, mass density

%% Be
% Z=4; % nuclear charge number
% Am=9.01218; % atomic weight
% M=Am*1e-3; % kg/mol, molar mass
% m=M/Na; % kg, mass of the atom

% Rho=1.85e3; % kg/m^3, mass density

%%
n=Rho/m; % m^-3, number density
aWS=(3./(4*pi.*n)).^(1/3); % m, Wigner-Seitz radius 
%% Temperature
T1=0.1:0.1:0.9;
T2=1:0.2:9;
T3=10:10:10000;
Te=[T1,T2,T3]; % eV, temperature range
Te=Te'; % eV, temperature
T=Te*elect/kb; % K, temperature
[NT,~]=size(T); % size of the data of T
%% Z_M, average ionization
%---- Thomas Fermi model -------
Z_M_TF = Z_TF(Z,Am,n,T); % TF model
Z_M=Z_M_TF; % average ionization
LD=sqrt(epsilon0*kb.*T./(Z_M*n)/elect^2./(Z_M+1)); % m, Debye length D,
mi=m-Z_M*me; % kg, ion mass when ionization equals to Z_M 
ne=Z_M*n; % m^-3, number density of electrons

%% Fermi energy
Ef=(h/2/pi)^2/2/me.*(3*pi^2.*Z*Rho/m).^(2/3); % J, Fermi Energy
Tf=Ef/kb; % K, Fermi temperature
Tfe=Tf/elect*kb; % eV, Fermi temperature

%% Calculate the intigrating term: I(T)
tic; % Timer start
I_T_ii=zeros(NT,1);
for j=1:NT
   I_T_ii(j,1)=Inti2ii(T(j),Z_M(j,1),LD(j,1)); % m^2, integral term

   % --------------- Timer ----------------
   percents=j/NT*100;
   clc;
   timepass=toc;
   fprintf('Calculating the inigrating term I(T): %.1f%% have been done.\n',percents);
   run_speed=j/timepass;
   fprintf('About %.0f min %.1f s remain.\n',floor((NT-j)/run_speed/60),mod((NT-j)/run_speed,60));
   % --------------- Timer ----------------
end

%% The viscosity
d=LD; % m, cut-off distance
%  ---- Fomular 7 ---- 
Eta_ii=1/pi./(d.^4).*sqrt(3*mi*kb.*T).*I_T_ii; % Pa s, viscosity
%  ---- Fomular 7 ---- 

Eta=Eta_ii*1000; % mPa s, viscosity

%% cut off data
% remove data at low temperatures
Nlower=1;
for icount=1:NT    
    if LD(icount)<0.1*aWS % cut-off condition
        Eta(icount)=nan;
        Nlower=Nlower+1; % counts
    end
end

% remove data at high temperatures
Nupper=NT;
for icount=1:NT
    if LD(icount)>0.5*aWS % cut-off condition
        Eta(icount)=nan;
        Nupper=Nupper-1;
    end
end
Tlower=Te(Nlower); % eV, lower limit of temperature
Etadown=Eta(Nlower); % mPa s, lower limit of viscosity
Tupper=Te(Nupper); % eV, upper limit of temperature
Etaupper=Eta(Nupper); % mPa s, upper limit of viscosity
Results_cutoff=[Tlower,Etadown,Tupper,Etaupper]; % eV, cut-off temperature and corresponding viscosity 

%% Figure
figure;
loglog(Te,Eta,'b');
xlabel('T (eV)');
ylabel('Eta (mPa s)');
legend('RWSP-VM');

%%
timepass=toc;
fprintf('Completed!\nTotal time cost: %.2f second(s)\n',timepass);




