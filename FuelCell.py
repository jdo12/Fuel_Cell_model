# -*- coding: utf-8 -*-
# !/usr/bin/env python3
# PEM FuelCell model v0.1.0
"""

"""
from math import log, exp
from numpy import zeros, convolve
from scipy.optimize import fsolve
from helpers import fit_pdat, vifc_calc
from matplotlib.pyplot import plot, title, show

# Model parameters to be included in config.ini
T = 70.0                                            # (°C) FuelCell operating temperature
Pe = 101000                                         # (pa)
Pe_out = 50000000                                   # (Pa) state of maximum charge in the tank
F = 96485.34                                        # (C/mol) Faraday's constant
A = 200                                             # (cm^2) area of a cell
ne = 2
charge_lvl_min = 19                                 # (%) minimum state of charge in the hydrogen tank
charge_lvl_max = 95                                 # (%) maximum state of charge in the hydrogen tank
soc_lvl_init = 95                                   # (%) initial state of charge
DG_25 = 237100.0
DG_200 = 220370.0
DH = 286000.0
R = 8.31447                                         # (J/mol-K) universal constant of gases
Ne = 50                                             # Number of Electrolyzers
V_tank = 0.3                                        # (m^3) volume of the tank
Nt = 8                                              # number of tanks
Nc = 20                                             # number of cells
Nfc = 50                                            # number of fuel cells
# V vs i characteristic curve
I_initial = 0
I_final = 200
I_step = 1
alpha = 0.5
a = 3e-5                                            # (V)
b = 8e-3*0.001                                      # (A)
i_n = 2.0*0.001                                     # (A)
i_0 = 0.067*0.001                                   # (A)
i_L = 900.0*0.001                                   # (A)
R_tot = 30e-6*1000                                  # (Ohm)
# Mass transport dynamics considered as time delay
Patm = 101000                                       #(Pa)
Tau_e = 80                                          # sec overall flow delay
Lambda = 0.00333                                    # (ohm) constant
pH2 = 0.5                                           # effective partial pressure of H2
pO2 = 0.8                                           # effective partial pressure of O2
x0_1 = 0.7
x0_2 = 80
power_data = 'pdata.csv'                            # Power curve time-series CSV data

# Compute starting parameters for model init
T += 273.15
min_charge = charge_lvl_min*1e-2*Pe_out             # min. state of charge in the hydrogen tank
max_charge = charge_lvl_max*1e-2*Pe_out             # max. state of charge in the hydrogen tank
soc_i = soc_lvl_init*1e-2*Pe_out                    # SOC initial
n_i = soc_i*V_tank/R/(25+273.15)                    # initial number of H2 moles in the tank
DG = DG_200+(200-T+273.15)/175*(DG_25-DG_200)
V_rev = DG/ne/F
Vtn = DH/ne/F                                       # thermo-neutral voltage

nn = int((I_final-I_initial)/I_step)

b1 = b/A                                            # (A)
i_n1 = i_n*A                                        # (A)
i_01 = i_0*A                                        # (A)
R_tot1 = R_tot/A                                    # (Ohm)

# Fit the power curve data
p, timespan = fit_pdat(power_data)

Pl, Pr, P_tot, Ir, P, ne, Qh2_V, Qh2_m,\
    m_dotH2, P_tank, Tout, Vi, Vr, = [zeros((timespan-1, 1), dtype=float) for i in range(13)]
moles, Videal, Vreal, soc, Vact, Vohm, Vcon = [zeros((timespan, 1), dtype=float) for i in range(7)]
id0, conv, expo, Ed_cell, E_cell, V_graph, Ii = [zeros((timespan-1, 1), dtype=float) for i in range(7)]

P_tank[0] = soc_i
moles[0][0] = n_i*Nt

for i in range(timespan-1):
    # Power profile (kW)
    Pl[i] = sum([p[j]*(i+1)**(21-j) for j in range(22)])
    Pr[i] = 1e3*Pl[i]/Nfc/Nc                        # (W) requested power for one cell

    Vi[i], Ir[i] = fsolve(vifc_calc, [x0_1, x0_2], args=(Pr[i], V_rev, T, R, alpha, F, i_n1, i_01, R_tot1, a, b1, pH2, pO2))

    Vact[i] = (R*T/alpha/F)*log((Ir[i]+i_n1)/i_01)
    Vohm[i] = (Ir[i]+i_n1)*R_tot1
    Vcon[i] = a*exp(b1*Ir[i])
    V_graph[i] = V_rev - Vact[i] - Vohm[i] + Vcon[i]
    Ii[i] = Ir[i]*1e3/A
    id0[i] = Ir[i]/A
    expo[i] = exp(i/Tau_e)
    conv[i] = convolve(id0[i], expo[i])
    Ed_cell[i] = Lambda*(id0[i]-conv[i])
    E_cell[i] = V_rev+R*T/2/F*log(pH2*pO2**0.5)-Ed_cell[i]
    Vr[i] = E_cell[i]-Vact[i]-Vohm[i]+Vcon[i]
    if Vr[i] > Vi[i]:
        Vr[i] = Vi[i]
    Videal[i] = Vi[i]
    Vreal[i] = Vr[i]

    P[i] = Vr[i]*Ir[i]*Nfc*Nc/1000                  # (kW)
    ne[i] = Vr[i]/1.48
    Qh2_m[i] = Nfc*Nc*Ir[i]/2/F                     # (mol/s)

    # Number of moles in time i in the tank
    moles[i+1] = moles[i]-Qh2_m[i]*1
    P_tank[i+1] = moles[i+1]*R*(T+273.15)/V_tank
    soc[i] = P_tank[i]

    if soc[i] < min_charge:
        print('Tank is fully discharged!!\n'
              'Discharge time:%d sec.' % i)

        soc_i = round((soc[i]/max_charge*100)[0], 1)
        print("State of Charge:%d%%" % soc_i)
        break

title('Ptank')
plot(P_tank[:i])
show()
title('Ir')
plot(Ir[:i])
show()
title('Vr')
plot(Vr[:i])
show()
title('moles')
plot(moles[:i])
show()
title('m_dotH2')
plot(m_dotH2[:i]*1000)                               # (g/s)
show()
title('Ptot')
plot(P_tot[:i])
show()
title('Pr')
plot(P[:i], 'r')
show()
title('Pr')
plot(P[:i], 'r')
show()
