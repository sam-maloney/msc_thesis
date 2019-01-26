gamma = 4/3; % R/Cv + 1
R = 5/7; % specific gas constant. depends on molWeight
Runiversal_foam = 8.31451e3;
Cv = R/(gamma-1); % 15/7
Cp = R + Cv; % 20/7
pin = 0.75;
rhoin = 1;
K = pin/rhoin^gamma; % equation of state coefficient
M1 = 5;
%Min = sqrt((2+(gamma-1)*M1^2)/(2*gamma*M1^2-gamma+1));
Min = sqrt(31/199);
cin = 1;
vin = -Min*cin;
v1 = vin*(gamma+1)*M1^2/(2+(gamma-1)*M1^2);
c1 = -v1/M1;
rho1 = rhoin*vin/v1;
p1 = rho1*c1^2/gamma;
Tin = cin^2/(gamma*R);
T1 = c1^2/(gamma*R);
Ber = 0.5*vin^2 + gamma/(gamma-1)*pin/rhoin + 0.935973200965766;
H = 1;
taac = H/abs(vin)/(1-Min);
omega0 = 4*pi/taac;
cout = sqrt(cin^2/0.75);
vout = vin*(cin/cout)^(2/(gamma-1));
Tout = cout^2/(gamma*R);
Mout = abs(vout/cout);
rhoout = rhoin*vin/vout;
pout = rhoout*cout^2/gamma;
dH = 0.1*H;
dPhi = (Mout^2/2 + 1/(gamma-1))*cout^2 - (Min^2/2 + 1/(gamma-1))*cin^2;

Lx = 4;
kx = 2*pi/Lx;
% muin = sqrt(1-kx^2*cin^2/omega0^2*(1-Min^2));
% muout = sqrt(1-kx^2*cout^2/omega0^2*(1-Mout^2));
% tQ = taac*(1+muin*Min)/(1+Min);
% Qnabla = (Mout+muout)/(1+muout*Mout)...
%          *exp(1i*omega0*tQ)/(muout*cin^2/cout^2 + muin*Mout/Min)...
%          *(1-cin^2/cout^2+kx^2*cin^2/omega0^2*(Min^2-Mout^2));
Qnabla_paper = 0.32; % approx as read from graph
Qnabla_sim = 0.338464957381226;

mu = sqrt(1-kx^2*cin^2/omega0^2*(1-Min^2));
kz_ = (omega0/cin)*(Min-mu)/(1-Min^2);
pert_mag = 1e-3*(1+mu*Min)/(1-Min^2);
dS_th = pert_mag*(2/Min)*(1-Min^2)/(1+gamma*Min^2)*(1-Min^2/M1^2)...
        *mu/(mu^2+2*mu*Min+M1^(-2));
dS_sim = 0.002940601;