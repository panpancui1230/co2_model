import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from calc import block

#parameters
lumen_protons_per_turnover = 0.000587
ATP_synthase_max_turnover = 200.0
Volts_per_charge = 0.047
perm_K = 150
n = 4.666
buffering_capacity = 0.03
k_KEA = 2500000
k_VCCN1 = 12
k_CLCE = 800000

def get_mutant_params(gtype):
    k_CLCE_mut = k_CLCE if gtype is None or 'clce2' not in gtype else 0
    k_KEA_mut = k_KEA if gtype is None or 'kea3' not in gtype else 0
    k_VCCN1_mut = k_VCCN1 if gtype is None or 'vccn1' not in gtype else 0
    return k_CLCE_mut, k_KEA_mut, k_VCCN1_mut


def model(y,t,k_CLCE_mut, k_KEA_mut, k_VCCN1_mut ):
    computer = block()

    pHlumen, Dy, pmf, Klumen, Kstroma, Cl_lumen, Cl_stroma, Hstroma, pHstroma=y

    #KEA3
    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    v_KEA = k_KEA_mut*(Hlumen*Kstroma -  Hstroma*Klumen)
   
    #V_K 可能要加调控
    K_deltaG=-0.06*np.log10(Kstroma/Klumen) + Dy
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2
    
    #VCCN1
    driving_force_Cl = 0.06* np.log10(Cl_stroma/Cl_lumen) + Dy
    v_VCCN1 = k_VCCN1_mut * computer.Cl_flux_relative(driving_force_Cl) * (Cl_stroma + Cl_lumen)/2

    #CLCE
    v_CLCE =  k_CLCE_mut*(driving_force_Cl*2+pmf)*(Cl_stroma + Cl_lumen)*(Hlumen+Hstroma)/4   

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


gtypes=[None,'clce2','kea3','vccn1']
labels=['WT','clce2','kea3','vccn1']
colors=['black','red','blue','green']

t = np.arange(0, 1200, 0.1)

plt.figure(figsize=(12, 8))

for gtype, label, color in zip(gtypes, labels, colors):
    k_CLCE_mut, k_KEA_mut, k_VCCN1_mut = get_mutant_params(gtype)
    sol = odeint(model, initial, t, args=(k_CLCE_mut, k_KEA_mut, k_VCCN1_mut))
    pmf = sol[:, 2]
    plt.plot(t, pmf, label=label, color=color)


# Plot settings
plt.xlabel('Time (s)')
plt.ylabel('pmf')
plt.title('Wild Type (WT) vs Mutants')
plt.legend()
plt.grid(True)
plt.show()
