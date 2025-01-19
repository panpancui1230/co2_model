import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

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
    def ATP_synthase_actvt(self, t):#based on gH+ data
        # x = t/T_ATP
        # actvt = 0.2 + 0.8*(x**4/(x**4 + 1))
        # return actvt
        actvt = 0.05 + 0.95 * np.exp(-0.008 * t)
        return actvt
    
    #合成
    def Vproton_pmf_actvt(self, pmf, actvt, ATP_synthase_max_turnover, n):# fraction of activity based on pmf, pmf_act is the half_max actvt pmf
        v_proton_active = 1 - (1 / (10 ** ((pmf - 0.162)*1.5/0.06) + 1))#reduced ATP synthase 2.7
        v_proton_inert = 1-(1 / (10 ** ((pmf - 0.204)*1.5/0.06) + 1))#oxidized ATP synthase 3.4
        
        v_active = actvt * v_proton_active * n * ATP_synthase_max_turnover
        v_inert = (1-actvt) * v_proton_inert * n * ATP_synthase_max_turnover
        
        v_proton_ATP = v_active + v_inert
        # v_proton_ATP = 0
        return (v_proton_ATP)
    

    #水解
    def V_ATP_proton_pmf_actvt(self, pmf, actvt, n):# fraction of activity based on pmf, pmf_act is the half_max actvt pmf
        v_ATP_proton_active = 1 - (1 / (10 ** ((pmf - 0.099)*1.5/0.06) + 1))#reduced ATP synthase 1.65
        v_ATP_proton_inert = 1-(1 / (10 ** ((pmf - 0.201)*1.5/0.06) + 1))#oxidized ATP synthase 3.35
        
        v_ATP_active = actvt * v_ATP_proton_active * n * 25
        v_ATP_inert = (1-actvt) * v_ATP_proton_inert * n * 25
        
        v_ATP_proton = v_ATP_active + v_ATP_inert
        # v_proton_ATP = 0
        return (v_ATP_proton)
        
    def V_H_dark(self, v_proton_ATP, v_ATP_proton, pmf, Hlumen, k_leak = 3*10**7):
    # def V_H_dark(self, v_proton_ATP, v_ATP_proton, pmf, Hlumen, k_leak = 0):
        V_H = v_proton_ATP - v_ATP_proton + pmf*k_leak*Hlumen
        # V_H = v_proton_ATP + pmf*k_leak*Hlumen
        return V_H

    #需要更新
    def Cl_flux_relative(self, v):
        # Cl_flux_v = 332*(v**3) + 30.8*(v**2) + 3.6*v
        Cl_flux_v = 29.044 * np.exp(648.585 * (v / 100)) -29.075
        return Cl_flux_v

def model(y,t):
    computer = block()

    pHlumen, Dy, pmf, Klumen, Kstroma, Cl_lumen, Cl_stroma, Hstroma, pHstroma=y

    #KEA3
    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    v_KEA = Kx.k_KEA*(Hlumen*Kstroma -  Hstroma*Klumen)
   
    #V_K
    K_deltaG= -0.06*np.log10(Kstroma/Klumen) + Dy
    # if Kstroma > 0 and Klumen > 0:
    #     K_deltaG = -0.06 * np.log10(Kstroma / Klumen) + Dy
    # else:
    #     K_deltaG = 0
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2 #修改
    # v_K_channel = 0 * K_deltaG*(Klumen+Kstroma)/2
    
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

    # dATP 
    activity = computer.ATP_synthase_actvt(t)
    d_protons_to_ATP = computer.Vproton_pmf_actvt(pmf, activity, ATP_synthase_max_turnover, n)
    d_ATP_to_protons = computer.V_ATP_proton_pmf_actvt(pmf, activity, n)
  
    #dpHlumen 
    d_H_ATP_or_passive = computer.V_H_dark( d_protons_to_ATP, d_ATP_to_protons, pmf, Hlumen,k_leak = 3*10**7)    
    net_protons_in =  - d_H_ATP_or_passive
    dHin = (net_protons_in - v_KEA - v_CLCE)*lumen_protons_per_turnover
    dpHlumen= -1*dHin / buffering_capacity 

    #dpHstroma 
    dHstroma = 0
    dpHstroma = -1*dHstroma / buffering_capacity

    #dDy 
    delta_charges= -d_H_ATP_or_passive - v_K_channel - v_VCCN1-3*v_CLCE 
    dDy=delta_charges*Volts_per_charge

    #dpmf #pHlumen正负号
    dpmf= -0.06* dpHlumen + dDy
    # dpmf= 0.06* dpHlumen + dDy

    return [dpHlumen, dDy, dpmf, dKlumen, dKstroma, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]

def sim_a_gtype(gtype='WT'):
    Kx.k_KEA = 2500000
    Kx.k_VCCN1 = 12
    Kx.k_CLCE = 800000

    if 'kea3' in gtype:
        Kx.k_KEA =0
    if 'vccn1' in gtype:
        Kx.k_VCCN1 =0
    if 'clce2' in gtype:
        Kx.k_CLCE =0

dpHlumen_initial = 5.994956
dDy_initial = 0.010221
dpmf_initial = 0.118524
dKlumen_initial = 0.086659
dKstrom_initial = 0.1
dCl_lumen_initial = 0.053123
dCl_stroma_initial = 0.038688
dHstroma_initial = 0.0
dpHstroma_initial = 7.8
initial=[dpHlumen_initial, dDy_initial, dpmf_initial, dKlumen_initial, dKstrom_initial, 
         dCl_lumen_initial, dCl_stroma_initial, dHstroma_initial, dpHstroma_initial]

t = np.arange(0,7200,0.1)
gtypes = ['WT', 'kea3', 'vccn1']
colors = ['black', 'blue', 'red']
variables = ['pHlumen', 'Dy', 'pmf', 'Klumen', 'Kstroma', 'Cl_lumen', 'Cl_stroma', 'Hstroma', 'pHstroma']

exclude_vars = ['Kstroma', 'Hstroma', 'pHstroma']  # Variables to exclude from plotting

# Determine the indices of variables to plot
plot_indices = [i for i, var in enumerate(variables) if var not in exclude_vars]
plot_vars = [variables[i] for i in plot_indices]

# Adjust font sizes globally
plt.rcParams.update({
    'font.size': 14,  # Set the base font size for all text
    'axes.titlesize': 20,  # Title size
    'axes.labelsize': 16,  # Label size
    'xtick.labelsize': 12,  # X-axis tick size
    'ytick.labelsize': 12,  # Y-axis tick size
    'legend.fontsize': 12  # Legend font size
})

fig, axes = plt.subplots(2, 3, figsize=(14, 8))  # Adjusted for 6 subplots
axes = axes.flatten()

for idx in range(len(gtypes)):
    gtype = gtypes[idx]
    color = colors[idx]
    
    sim_a_gtype(gtype)
    sol = odeint(model, initial, t)
    
    for i, var_idx in enumerate(plot_indices):
        axes[i].plot(t, sol[:, var_idx], label=gtype, color=color, alpha=0.75)

for i, var in enumerate(plot_vars):
    # axes[i].set_title(var)
    axes[i].set_xlabel('Time/s')
    axes[i].set_ylabel(var)
    # axes[i].legend(loc='upper right')
    axes[i].legend(
        loc='upper center', 
        bbox_to_anchor=(0.5, 1.1),  # Adjust position above the plot
        ncol=3,                    # Number of columns for the legend
        frameon=False,             # No border around the legend
        fontsize=10                # Smaller font size for clarity
    )
    axes[i].grid(False)

plt.tight_layout()
plt.show()





# gtypes = ['WT', 'kea3', 'vccn1']
# # gtypes = ['kea3']
# colors = ['black', 'blue', 'red']
# # colors = ['blue']
# variables = ['pHlumen', 'Dy', 'pmf', 'Klumen', 'Kstroma', 'Cl_lumen', 'Cl_stroma', 'Hstroma', 'pHstroma']

# results = {}
# for gtype in gtypes:
#     sim_a_gtype(gtype)
#     sol = odeint(model, initial, t)
#     results[gtype] = sol

# delta_K = {gtype: results[gtype][:,3] - 0.1 for gtype in gtypes}
# delta_Cl = {gtype: results[gtype][:,5] - results[gtype][:,6]
#             for gtype in gtypes}
# # lumen_total = {gtype: results[gtype][:,3] + results[gtype][:,5]
# #                for gtype in gtypes}
# # stroma_total = {gtype: results[gtype][:,6] + 0.1 for gtype in gtypes}
# # delta_total = {gtype: lumen_total[gtype] - stroma_total[gtype]
# #                for gtype in gtypes}

# # for gtype in gtypes:
# #     print(f"Genotype: {gtype}")
# #     print(f"Lumen Total (First 10): {lumen_total[gtype][:10]}")
# #     print(f"Stroma Total (First 10): {stroma_total[gtype][:10]}")
# #     print(f"Delta Total (First 10): {delta_total[gtype][:10]}")

# plt.figure(figsize=(14,6))

# # for idx, gtype in enumerate(gtypes):
# #     plt.plot(t,delta_total[gtype], label = gtype, color = colors[idx], alpha = 0.75)
    
# # plt.xlabel('Time(s)', fontsize = 16)
# # plt.ylabel('delta_total (lumen_total-stroma_total)',  fontsize = 16)
# # plt.legend(
# #     loc = 'upper center',
# #     bbox_to_anchor = (0.5,1.1),
# #     ncol = 3,
# #     frameon = False,
# #     fontsize = 18
# # )
# plt.subplot(1,2,1)
# for idx, gtype in enumerate(gtypes):
#     plt.plot(t, delta_K[gtype], label = gtype, color = colors[idx], alpha = 0.75)
# plt.xlabel('Time(s)', fontsize=16)
# plt.ylabel('delta_K (Klumen-0.1)', fontsize=16)
# # plt.title
# plt.legend(
#     loc = "upper center",
#     bbox_to_anchor = (0.5,1.1),
#     ncol = 3,
#     frameon = False,
#     fontsize = 12
# )
# plt.grid(False)

# plt.subplot(1,2,2)
# for idx, gtype in enumerate(gtypes):
#     plt.plot(t, delta_Cl[gtype], label = gtype, color = colors[idx], alpha = 0.75)
# plt.xlabel('Time(s)', fontsize=16)
# plt.ylabel('delta_Cl (Cl_lumen-Cl_stroma)', fontsize=16)
# # plt.title
# plt.legend(
#     loc = "upper center",
#     bbox_to_anchor = (0.5,1.1),
#     ncol = 3,
#     frameon = False,
#     fontsize = 12
# )
# plt.grid(False)

# plt.show()


