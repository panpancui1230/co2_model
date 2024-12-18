import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import solve_ivp
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import importlib as im
from matplotlib import cm
import copy
import pandas as pd
from scipy import integrate
from scipy import signal
from IPython.core.display import display, HTML
import csv
import warnings

from painter import Plotting
from utils import standard_constants
from utils import standard_initial_states
from utils import sim_states
from utils import sim_constants
from calc import block
# from sun_sim import sunshine

def f(t, y, pKreg, max_b6f, lumen_protons_per_turnover, ATP_synthase_max_turnover, Volts_per_charge, perm_K,
      n, Em7_PQH2, Em7_PC, Em_Fd, buffering_capacity, k_Fd_to_NADP, k_CBC, k_KEA, k_VCCN1, k_CLCE, k_NDH):
    computer = block()
    y = Dy, pHstroma, pHlumen, 
    pmf=Dy + 0.06*(pHstroma-pHlumen)
    b6f_content=0.433 

    v_b6f=computer.calc_v_b6f(max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf)
    
    v_NDH = computer.calc_v_NDH(Em_Fd, Em7_PQH2, pHstroma, pmf, k_NDH, Fd_red, Fd_ox, PQ, PQH2)
    d_Hlumen_NDH = v_NDH*2 
    d_charge_NDH = d_Hlumen_NDH 

    PGR_vmax = 0
    v_PGR = computer.calc_v_PGR(PGR_vmax, Fd_red, PQ, PQH2)

    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    activity = computer.ATP_synthase_actvt(t, T_ATP)

    d_protons_to_ATP = computer.Vproton_pmf_actvt(pmf, activity, ATP_synthase_max_turnover, n)
    d_H_ATP_or_passive = computer.V_H_light(light_per_L, d_protons_to_ATP, pmf, Hlumen)                              

        
    d_ATP_made=d_protons_to_ATP/n 

    NADPH_CBC = k_CBC*(1.0-np.exp(-t/900))*(np.log(NADPH_pool/NADP_pool)-np.log(1.25))/(np.log(3.5/1.25))#calc_CBC_NADPH(k_CBC, t, d_ATP_made)
   
    dNADPH_pool=0.5 * k_Fd_to_NADP*NADP_pool*Fd_red - NADPH_CBC
    dNADP_pool=-1*dNADPH_pool
    
    d_ATP_consumed = d_ATP_made

    d_protons_from_b6f = v_b6f*2 

    charges_from_b6f = v_b6f

    net_protons_in = d_protons_from_PSII + d_protons_from_b6f + d_Hlumen_NDH - d_H_ATP_or_passive
   
    
    f_actvt = computer.KEA_reg(pHlumen, QAm)
    v_KEA = k_KEA*(Hlumen*Kstroma -  Hstroma*Klumen)*f_actvt#/(10**(2*(pHlumen-6.5))+1)
    
    #Pck = 1/(1+np.exp(39.5*0.66*(-0.003-Dy)))#probability of v_K_channel open
    K_deltaG=-0.06*np.log10(Kstroma/Klumen) + Dy
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2
    
    net_Klumen =  v_KEA - v_K_channel        
    dKlumen = net_Klumen*lumen_protons_per_turnover   
    dKstroma=0

    driving_force_Cl = 0.06* np.log10(Cl_stroma/Cl_lumen) + Dy
    v_VCCN1 = k_VCCN1 * computer.Cl_flux_relative(driving_force_Cl) * (Cl_stroma + Cl_lumen)/2#Cl_flux_relative()

    v_CLCE =  k_CLCE*(driving_force_Cl*2+pmf)*(Cl_stroma + Cl_lumen)*(Hlumen+Hstroma)/4    
    net_Cl_lumen_in = v_VCCN1 + 2*v_CLCE
    dCl_lumen = net_Cl_lumen_in * lumen_protons_per_turnover
    dCl_stroma = -0.1*dCl_lumen
    

    dHin = (net_protons_in - v_KEA - v_CLCE)*lumen_protons_per_turnover
    dpHlumen= -1*dHin / buffering_capacity 

    dHstroma = 0
    
    dpHstroma = -1*dHstroma / buffering_capacity
    delta_charges=charges_from_PSII+PSI_charge_separations + charges_from_b6f \
                    + d_charge_NDH - v_K_channel - d_H_ATP_or_passive - v_VCCN1-3*v_CLCE
    
    dDy=delta_charges*Volts_per_charge
    dpmf= 0.06* dpHlumen + dDy

    dATP_pool= d_ATP_made - d_ATP_consumed
    dADP_pool= - dATP_pool
    ddeltaGatp = 0
    return [ dHin, dpHlumen, dDy, dpmf, ddeltaGatp, dKlumen, dKstroma, 
            d_ATP_made, dATP_pool, dADP_pool, dNADPH_pool,dNADP_pool, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]

