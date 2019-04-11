% PEM fuel cell

% fuel cell parameters
T=70+273.15; % (°C)operating temperature
Pe=101000; % (pa)
Pe_out = 50000000; % (Pa) state of maximun charge in the tank
F = 96485.34; % (C/mol) Faraday's constant
ne = 2; % 
DG_25 = 237100;
DG_200 = 220370;
DG = DG_200+(200-T+273.15)/175*(DG_25-DG_200);
Vrev = DG/ne/F;
DH = 286000;
A = 200; % (cm^2) area of a cell
Vtn = DH/ne/F; % thermoneutral voltage
min_charge = 0.19*Pe_out; % minumun state of charge in the hydrogen tank
max_charge = 0.95*Pe_out; % maximun state of charge in the hydrogen tank
R = 8.31445; % (J/mol-K) universal constant of gases
Vtank = 0.3; % (m^3) volume of the tank

Nc = 20; %(number of cells)
Nfc = 50; % number of fuel cells
Nt = 8; % number of tanks
%Nt = 5000; % number of tanks % for aging test only

tvector = importdata('time.txt');
Pvector = importdata ('power.txt');
X=tvector;
Y=Pvector;
p=polyfit(X,Y,21);
timespan = X(length(X));
%timespan = 5040000; % for aging test only

Pl=zeros(timespan-1,1);
%Pr=zeros(timespan-1,1);
Ptot=zeros(timespan-1,1);
Ir=zeros(timespan-1,1);
P=zeros(timespan-1,1);
V=zeros(timespan-1,1);
n_e=zeros(timespan-1,1);
Qh2_V=zeros(timespan-1,1);
Qh2_m=zeros(timespan-1,1);
m_dotH2=zeros(timespan-1,1);
Ptank=zeros(timespan-1,1);
Tout=zeros(timespan-1,1);
Vi=zeros(timespan-1,1);
Vr=zeros(timespan-1,1);
moles=zeros(timespan,1);
Videal=zeros(timespan,1);
Vreal=zeros(timespan,1);
Soc=zeros(timespan,1);
Vact=zeros(timespan,1);
Vohm=zeros(timespan,1);
Vcon=zeros(timespan,1);
t=zeros(timespan,1);

% V vs i characteristic curve
Iinitial = 0;
Ifinal = 200; % (A)
Istep = 1;
nn = (Ifinal-Iinitial)/Istep;

alpha = 0.5;
a = 3e-5; % (V)
b = 8e-3*0.001/A; % (A)
i_n = 2.0*0.001*A; % (A)
i_0 = 0.067*0.001*A; % (A)
i_L = 900.0*0.001*A; %(A)
R_tot = 30e-6*1000/A; % (Ohm)

% Mass transport dynamics considered as time delay
Patm = 101000; %(Pa)
Tau_e = 80; % sec overall flow dealy
Lambda = 0.00333; %(ohm) constant
pH2 = 0.5; % effective partial pressure of H2
pO2 = 0.8; % effective partial pressure of O2
id =zeros(nn+1,1);
Conv = zeros(nn+1,1);
expo = zeros(nn+1,1);
Ed_cell = zeros(nn+1,1);
E_cell = zeros(nn+1,1);

ii=[Iinitial:Istep:Ifinal]'; % A
Vgraph1=zeros(nn+1,1);
Vgraph=zeros(nn+1,1);
Ii=zeros(nn+1,1);

% for i = 1:nn+1
%     Vact(i) = (R*T/alpha/F)*log((ii(i)+i_n)/i_0);
%     Vohm(i) = (ii(i)+i_n)*R_tot;
%     %Vcon = -(R*T/2/F)*log(1-(ii(i)+i_n)/i_L);
%     Vcon(i) = a*exp(b*ii(i));
%     Vgraph1(i) = Vrev - Vact(i) - Vohm(i) + Vcon(i);
%     Ii(i) = ii(i)*1000/A;
%     id(i) = ii(i)/A;
%     expo(i) = exp(i/Tau_e);
%     Conv(i)=conv(id(i),expo(i));
%     Ed_cell(i) = Lambda*(id(i)-Conv(i));
%     E_cell(i) = Vrev + R*T/2/F*log(pH2*(pO2)^0.5)-Ed_cell(i);
%     Vgraph(i) = E_cell(i)- Vact(i) - Vohm(i) + Vcon(i);
% end
% plot(Ii,Vgraph,Ii,Vgraph1,'r')
%3234392 % current iteration for metrics

Ptank(1)=0.95*Pe_out; % initial state of charge 
Ptank_age(1)=0.95*Pe_out; % initial state of charge %%%%%%%%% aging
n_i = Ptank(1)*0.3/8.3145/(25+273.15); % initial number of H2 moles in the tank 
n_i_age = Ptank_age(1)*0.3/8.3145/(25+273.15); % initial number of H2 moles in the tank %%%%%%%%% aging
moles(1)=n_i*Nt; %%%%%%%%% changed
moles_age(1)=n_i_age*Nt; %%%%%%%%% aging
Lka_H2(1)=0; %%%%%%%% leakage
for i = 1:timespan-1
t(i)=i;

%Power profile (kW)
Pl(i) = (p(1)*t(i)^21+p(2)*t(i)^20+p(3)*t(i)^19+p(4)*t(i)^18+p(5)*t(i)^17+p(6)*t(i)^16+p(7)*t(i)^15+p(8)*t(i)^14+p(9)*t(i)^13+p(10)*t(i)^12+p(11)*t(i)^11+p(12)*t(i)^10+p(13)*t(i)^9+p(14)*t(i)^8+p(15)*t(i)^7+p(16)*t(i)^6+p(17)*t(i)^5+p(18)*t(i)^4+p(19)*t(i)^3+p(20)*t(i)^2+p(21)*t(i)+p(22));
Pr(i) = 1000*Pl(i)/Nfc/Nc; % (W) requested power for one cell

Pri=Pr(i);
 % Pr(i) = 100; % for aging test only
 % Pri = 100; % for aging test only

x0=[0.7 80];
x = fsolve(@(x)VIFC(x,Pri,Vrev,T,R,alpha,F,i_n,i_0,R_tot,a,b,pH2,pO2),x0);
Vi(i)=x(1);
Ir(i)=x(2);
Vact(i) = (R*T/alpha/F)*log((Ir(i)+i_n)/i_0);
Vohm(i) = (Ir(i)+i_n)*R_tot;
Vcon(i) = a*exp(b*Ir(i));
Ii(i) = Ir(i)*1000/A;
id(i) = Ir(i)/A;
expo(i) = exp(i/Tau_e);
Conv(i)=conv(id(i),expo(i));
Ed_cell(i) = Lambda*(id(i)-Conv(i));
E_cell(i) = Vrev + R*T/2/F*log(pH2*(pO2)^0.5)-Ed_cell(i);
Vr(i) = E_cell(i)- Vact(i) - Vohm(i) + Vcon(i);
 if Vr(i)> Vi(i)
    Vr(i)= Vi(i);
 end
 Videal(i) = Vi(i);
 Vreal(i) = Vr(i);
 Vr_age(i) = Vr(i)-(5.28489722540435e-15)*t(i)^2-(5.36476231386010e-08)*t(i); %%%%%%%%% aging
 Ir_age(i) = Pr(i)./Vr_age(i);
 
P(i)= Vr(i)*Ir(i)*Nfc*Nc/1000; % (kW)
P_age(i)= Vr_age(i)*Ir_age(i)*Nfc*Nc/1000; % (kW) %%%%%%%%% aging
n_e(i)=Vr(i)/1.48;
n_e_age(i)=Vr_age(i)/1.48; %%%%%%%%% aging

% Flow of H2 consumed
Qh2_m(i)= Nfc*Nc*Ir(i)/2/F; % (mol/s)
Qh2_m_age(i)= Nfc*Nc*Ir_age(i)/2/F; % (mol/s)%%%%%%%%% aging
 m_dotH2(i)=Qh2_m(i)*2/1000; % (kg/s)
 m_dotH2_age(i)=Qh2_m_age(i)*2/1000; % (kg/s)  %%%%%%%%%% aging

% Storage tank
moles(i+1) = moles(i) - Qh2_m(i)*1; % number of moles in time i in the tank(s)
Ptank(i+1)= moles(i+1)/Nt*R*(25+273.15)/Vtank;


  % hydrogen leakage Hydrogen Tank 525bar – 300L (MAHYTEC)- 403 Stainless steel 
  At = 47*Nt; % m^2 internal surface area of the tank(s) %%%%%%%% leakage
  l = 2.54/1000; % m tank wall thickness (2.54 mm) %%%%%%%% leakage
  Phi_0 = 5.9e-5; % mol s^-1 m^-1 MPa^-0.5 %%%%%%%% leakage
  E_phi = 42.7; % KJ/mol %%%%%%%% leakage
  b = 1.55e-5; % m^3/mol co-volume factor %%%%%%%% leakage
  Phi = Phi_0*exp(-E_phi/(R/1000)/(25+273.15)); % mol s^-1 m^-1 MPa^-0.5 %%%%%%%% leakage
  f(i)= Ptank(i)/1e6*exp(Ptank(i)/1e6*b/(R/1000)/(25+273.15)); % fugacity %%%%%%%% leakage
  J(i) = Phi/l*2*(f(i))^0.5; % g m^-2 s^-1 %%%%%%%% leakage
  Lk_H2(i) = J(i)*At; % g/s %%%%%%%% leakage
  Lka_H2(i+1)= Lka_H2(i) + J(i)*At*1; % total hydrogen lost in time t (g) %%%%%%%% leakage
  
  moles_age(i+1) = moles_age(i) - Qh2_m_age(i)*1-Lka_H2(i+1)/2*1; % number of moles in time i in the tank %%%%%%%%%% aging and leakage
  Ptank_age(i+1)= moles_age(i+1)/Nt*R*(25+273.15)/Vtank;%%%%%%%%% aging
  Soc_age(i) = Ptank_age(i)/max_charge; %%%%%%%%% aging
  
  % Discharging efficiency
  LHV_H2 = 120000; % (J/g)
  Eta_ds(i) = P(i)*1000/LHV_H2/(m_dotH2(i)*1000); %(W/W)

  if Ptank(i+1)< min_charge 
      disp('Tank is fully discharged charge:')
      disp(['discharging_time: ', num2str(i),' seconds']) 
      disp(['State of charge ',num2str(Soc(i)*100),'%']) 
  break
  end
end

%   plot(1:i,Ptank(1:i));
%   plot(1:i,Ir(1:i));
%   plot(1:i,Vr(1:i));
%   plot(1:i,moles(1:i));
%   plot(m_dotH2(1:i)*1000); % (g/s)
%   plot(1:i,ne(1:i));
%   plotyy(1:i,Ptot(1:i),1:i,Pr(1:i))

Results = cat(2,Ptank(1:i),Ir(1:i),Vi(1:i),moles(1:i),Qh2_m(1:i),n_e(1:i));