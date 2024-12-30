import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# from calc import block

lumen_protons_per_turnover = 0.000587
ATP_synthase_max_turnover = 200.0
Volts_per_charge = 0.047
perm_K = 150
n = 4.666
buffering_capacity = 0.03
class Kx:
    k_KEA = 2500000
    k_VCCN1 = 12
    k_CLCE = 800000

class block:
    def ATP_synthase_actvt(self, t, T_ATP):#based on gH+ data
        x = t/T_ATP
        actvt = 0.2 + 0.8*(x**4/(x**4 + 1))
        return actvt
    
    def Vproton_pmf_actvt(self, pmf, actvt, ATP_synthase_max_turnover, n):# fraction of activity based on pmf, pmf_act is the half_max actvt pmf
        v_proton_active = 1 - (1 / (10 ** ((pmf - 0.132)*1.5/0.06) + 1))#reduced ATP synthase
        v_proton_inert = 1-(1 / (10 ** ((pmf - 0.204)*1.5/0.06) + 1))#oxidized ATP synthase
        
        v_active = actvt * v_proton_active * n * ATP_synthase_max_turnover
        v_inert = (1-actvt) * v_proton_inert * n * ATP_synthase_max_turnover
        
        v_proton_ATP = v_active + v_inert
        return (v_proton_ATP)

    def V_H_dark(self, v_proton_ATP, pmf, Hlumen, k_leak = 3*10**7):
        V_H = -v_proton_ATP + pmf*k_leak*Hlumen
        return V_H

    def Cl_flux_relative(self, v):
        Cl_flux_v = 332*(v**3) + 30.8*(v**2) + 3.6*v
        #relative to Cl flux thru VCCN1. when driving force is 0.1 Volt,
        #Cl_flux_v is 1. empirical equation was obtained from
        # Herdean et al. 2016 DOI: 10.1038/ncomms11654
        return Cl_flux_v

def model(y,t):
    computer = block()

    pHlumen, Dy, pmf, Klumen, Kstroma, Cl_lumen, Cl_stroma, Hstroma, pHstroma=y

    #KEA3
    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    v_KEA = Kx.k_KEA*(Hlumen*Kstroma -  Hstroma*Klumen)
   
    #V_K 可能要加调控
    K_deltaG=-0.06*np.log10(Kstroma/Klumen) + Dy
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2
    
    #VCCN1
    driving_force_Cl = 0.06* np.log10(Cl_stroma/Cl_lumen) + Dy
    v_VCCN1 = Kx.k_VCCN1 * computer.Cl_flux_relative(driving_force_Cl) * (Cl_stroma + Cl_lumen)/2

    #CLCE
    v_CLCE =  Kx.k_CLCE*(driving_force_Cl*2+pmf)*(Cl_stroma + Cl_lumen)*(Hlumen+Hstroma)/4   

    #dK+/dt
    # net_Klumen =  v_KEA - v_K_channel        
    dKlumen = (v_KEA - v_K_channel)*lumen_protons_per_turnover  
    dKstroma=0

    #dCl-/dt
    # net_Cl_lumen_in = v_VCCN1 + 2*v_CLCE
    dCl_lumen = (v_VCCN1 + 2*v_CLCE) * lumen_protons_per_turnover
    dCl_stroma = -0.1*dCl_lumen

    # dATP 水解
    d_protons_to_ATP = computer.Vproton_pmf_actvt(pmf, 0, ATP_synthase_max_turnover, n)

    #dpHlumen 质子浓度    
    d_H_ATP_or_passive = computer.V_H_dark(d_protons_to_ATP, pmf, Hlumen)    
    net_protons_in =  - d_H_ATP_or_passive
    dHin = (net_protons_in - v_KEA - v_CLCE)*lumen_protons_per_turnover
    dpHlumen= -1*dHin / buffering_capacity 

    #dpHstroma 质子浓度
    dHstroma = 0
    dpHstroma = -1*dHstroma / buffering_capacity

    #dDy 电位差
    delta_charges= - v_K_channel - v_VCCN1-3*v_CLCE 
    dDy=delta_charges*Volts_per_charge

    #dpmf 质子驱动力
    dpmf= 0.06* dpHlumen + dDy

    return [dpHlumen, dDy, dpmf, dKlumen, dKstroma, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]

def sim_a_gtype(gtype='WT'):
    Kx.k_KEA = 2500000
    Kx.k_VCCN1 = 12
    Kx.k_CLCE = 800000

    if gtype == 'kea3':
        Kx.k_KEA = 0
    if 'kea3' in gtype:
        Kx.k_KEA =0
    if 'vccn1' in gtype:
        Kx.k_VCCN1 =0

dpHlumen_initial = 6.5
dDy_initial = 0.056
dpmf_initial = 0.112
dKlumen_initial = 0.1
dKstrom_initial = 0.1
dCl_lumen_initial = 0.04
dCl_stroma_initial = 0.04
dHstroma_initial = 0.0
dpHstroma_initial = 7.8
initial=[dpHlumen_initial, dDy_initial, dpmf_initial, dKlumen_initial, dKstrom_initial, 
         dCl_lumen_initial, dCl_stroma_initial, dHstroma_initial, dpHstroma_initial]

t = np.arange(0,1200,0.1)

gtypes = ['WT', 'kea3', 'vccn1', 'clce2']
colors = ['black', 'blue', 'green', 'red']
variables = ['pHlumen', 'Dy', 'pmf', 'Klumen', 'Kstroma', 'Cl_lumen', 'Cl_stroma', 'Hstroma', 'pHstroma']

fig, axes = plt.subplots(3, 3, figsize=(18, 12))
axes = axes.flatten()

for idx in range(len(gtypes)):
    gtype = gtypes[idx]
    color = colors[idx]
    
    sim_a_gtype(gtype)
    sol = odeint(model, initial, t)
    
    for i, var in enumerate(variables):
        axes[i].plot(t, sol[:, i], label=gtype, color=color)

for i, var in enumerate(variables):
    axes[i].set_title(var)
    axes[i].set_xlabel('Time (s)')
    axes[i].set_ylabel(var)
    axes[i].legend(loc='upper right')
    axes[i].grid(True)

plt.tight_layout()
plt.show()