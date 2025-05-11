import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import odeint
from scipy.integrate import solve_ivp

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
        # actvt = 0.05 + 0.95 * np.exp(-0.008 * t)
        actvt = 0.05 + 0.95 * np.exp(-0.00154 * t)
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
        # v_ATP_active = actvt * v_ATP_proton_active * n * 250
        # v_ATP_inert = (1-actvt) * v_ATP_proton_inert * n * 250
        
        v_ATP_proton = v_ATP_active + v_ATP_inert
        # v_proton_ATP = 0
        return (v_ATP_proton)
        
    def V_H_dark(self, v_proton_ATP, v_ATP_proton, pmf, Hlumen, k_leak = 3*10**7):
    # def V_H_dark(self, v_proton_ATP, v_ATP_proton, pmf, Hlumen, k_leak = 0):
        V_H = v_proton_ATP - v_ATP_proton + pmf*k_leak*Hlumen
        # V_H = -v_proton_ATP + pmf*k_leak*Hlumen
        return V_H
    
    #需要更新
    def Cl_flux_relative(self, v):
        # Cl_flux_v = 332*(v**3) + 30.8*(v**2) + 3.6*v
        Cl_flux_v = 29.044 * np.exp(648.585 * (v / 100)) -29.075
        return Cl_flux_v
    #Hagino,T.et al.Cryo-EM structures of thylakoid-located voltage-dependent chloride channel VCCN1.Nature Communications, (2022) 13:2505
def model(t,y):
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
    # dKstroma= -0.1*dKlumen
    dKstroma= 0
    
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

# t = np.arange(0,7200,0.1)
t_span = (0,7200)
t_eval = np.linspace(0,7200,7200)

results = {}
for gtype in gtypes:
    sim_a_gtype(gtype)
    # sol = odeint(model, initial, t)
    sol = solve_ivp(model, t_span, initial, t_eval = t_eval, method = "BDF")
    results[gtype] = sol

delta_K = {gtype: results[gtype].y[3,:] - 0.1 for gtype in gtypes}
delta_Cl = {gtype: results[gtype].y[5,:] - results[gtype].y[6,:]
            for gtype in gtypes}
lumen_total = {gtype: results[gtype].y[3,:] + results[gtype].y[5,:]
               for gtype in gtypes}
stroma_total = {gtype: results[gtype].y[6,:] + 0.1 for gtype in gtypes}
delta_total = {gtype: lumen_total[gtype] - stroma_total[gtype]
               for gtype in gtypes}

def parameter_individual_plots():
    """
    Run simulations for three groups of parameter combinations and generate a separate plot for each configuration.
    Groups:
    1. k_VCCN1=12, ATP_hydrolysis_turnover=25, k_KEA3=[0, 2500000, 5000000, 12500000, 25000000]
    2. k_KEA3=2500000, ATP_hydrolysis_turnover=25, k_VCCN1=[0, 12, 24, 60, 120]
    3. k_VCCN1=12, k_KEA3=2500000, ATP_hydrolysis_turnover=[0, 25, 50, 125, 250]
    Each plot shows delta_K and delta_Cl for a single parameter configuration.
    """
    # Parameter groups
    param_groups = [
        {
            'name': 'k_KEA3',
            'values': [0, 2500000, 5000000, 12500000, 25000000],
            'fixed': {'k_VCCN1': 12, 'ATP_hydrolysis_turnover': 25},
            'group_id': 'group1'
        },
        {
            'name': 'k_VCCN1',
            'values': [0, 12, 24, 60, 120],
            'fixed': {'k_KEA3': 2500000, 'ATP_hydrolysis_turnover': 25},
            'group_id': 'group2'
        },
        {
            'name': 'ATP_hydrolysis_turnover',
            'values': [0, 25, 50, 125, 250],
            'fixed': {'k_VCCN1': 12, 'k_KEA3': 2500000},
            'group_id': 'group3'
        }
    ]

    # Run simulations and generate plots for each group
    for group in param_groups:
        param_name = group['name']
        param_values = group['values']
        fixed_params = group['fixed']
        group_id = group['group_id']
        
        for param_value in param_values:
            # Set fixed parameters
            Kx.k_KEA = fixed_params.get('k_KEA3', Kx.k_KEA)
            Kx.k_VCCN1 = fixed_params.get('k_VCCN1', Kx.k_VCCN1)
            ATP_hydrolysis_turnover = fixed_params.get('ATP_hydrolysis_turnover', 25)
            
            # Set variable parameter
            if param_name == 'k_KEA3':
                Kx.k_KEA = param_value
            elif param_name == 'k_VCCN1':
                Kx.k_VCCN1 = param_value
            elif param_name == 'ATP_hydrolysis_turnover':
                ATP_hydrolysis_turnover = param_value
            
            # Run simulation
            sol = solve_ivp(model, t_span, initial, t_eval=t_eval, method="BDF", 
                           args=(ATP_hydrolysis_turnover,))
            
            # Calculate delta_K and delta_Cl
            delta_K = sol.y[3, :] - 0.1
            delta_Cl = sol.y[5, :] - sol.y[6, :]
            
            # Create a new figure for this parameter configuration
            plt.figure(figsize=(14, 6))
            # Plot delta_K
            plt.subplot(1, 2, 1)
            plt.plot(t_eval, delta_K, color='black', alpha=0.75)
            plt.xlabel('Time (s)', fontsize=16)
            plt.ylabel('ΔK (Klumen - 0.1)', fontsize=16)
            plt.title(f'{param_name} = {param_value}', fontsize=16)
            plt.grid(False)
            
            # Plot delta_Cl
            plt.subplot(1, 2, 2)
            plt.plot(t_eval, delta_Cl, color='black', alpha=0.75)
            plt.xlabel('Time (s)', fontsize=16)
            plt.ylabel('ΔCl (Cl_lumen - Cl_stroma)', fontsize=16)
            plt.title(f'{param_name} = {param_value}', fontsize=16)
            plt.grid(False)
            
            plt.tight_layout()
            # Save with a unique filename
            plt.savefig(f'param_{group_id}_{param_name}_{param_value}.png')
            plt.close()

# Run the parameter individual plots
parameter_individual_plots()