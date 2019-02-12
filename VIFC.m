 function F = VIFC(x,Pri,Vrev,T,R,alpha,F,i_n,i_0,R_tot,a,b,pH2,pO2)
 
Vr = x(1);
Ir = x(2);
Vact = (R*T/alpha/F)*log((Ir+i_n)/i_0);
Vohm = (Ir+i_n)*R_tot;
Vcon = a*exp(b*Ir);

E_cell = Vrev + R*T/2/F*log(pH2*(pO2)^0.5);

fa = Vr - E_cell + Vact + Vohm - Vcon;
fb = Vr*Ir-Pri;

F=[fa;fb];
 end