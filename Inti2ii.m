%% vicinity collision
function I_T=Inti2ii(T,Z_M,LD)
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
qs=(Z_M*elect).^2/4/pi/epsilon0;

%% parameters of the hyperbola
% equation: x^2/a^2-y^2/b^2=1

r0=LD; % m, cut off distance

a=qs./(3*kb*T+2*qs/r0); % m, semi-major axis
db=1e-14; % m, integration step
bm=sqrt(r0^2-2*a*r0); % m, upper limit of the integral
b=db:db:bm; % m, semi-minor axis
b=b'; % transpose operation
c=sqrt(a.^2+b.^2); % m, focal length
%% coordinate when r=r0
x0=a./c.*(r0-a); % m 
y0=b./c.*sqrt((r0-a).^2-c.^2); % m

% tangent equation at (x0,y0): y=k1*x+k2
k1=x0.*b.^2./(y0.*a.^2); % slope

%% angle between velocity and x-axis
theta=atan(k1); % rad
%% intigral term
I_term=b.*(sin(pi-2*theta)).^2; % m
%% ingegral results
I_T=sum(I_term)*db; % m^2

end




