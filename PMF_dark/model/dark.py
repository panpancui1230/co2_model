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

# from painter import Plotting
from utils import standard_constants
from utils import standard_initial_states
from utils import sim_states
from utils import sim_constants
from calc import block
# from sun_sim import sunshine


species_labels = [
    'QA', # 0
    'QAm', #1 
    'PQ', #2
    'PQH2', #3
    'Hin', #4
    'pHlumen', #5
    'Dy', #6
    'pmf', #7
    'DeltaGatp', #8
    'Klumen', #9
    'Kstroma', #10
    'ATP_made', #11
    'PC_ox', #12
    'PC_red', #13
    'P700_ox', #14
    'P700_red', #15
    'Z_array', #16
    'V_array', #17
    'NPQ_array', #18
    'singletO2_array', #19
    'Phi2_array', #20
    'LEF_array', #21
    'Fd_ox',
    'Fd_red',
    'ATP_pool',
    'ADP_pool',
    'NADPH_pool',
    'NADP_pool',
    'Cl_lumen',
    'Cl_stroma',
    'Hstroma',
    'pHstroma'
    ]


def f(t, y, lumen_protons_per_turnover, ATP_synthase_max_turnover, Volts_per_charge, perm_K,
      n, buffering_capacity, k_KEA, k_VCCN1, k_CLCE):
    computer = block()

    pHlumen, Dy, pmf, Klumen, Kstroma, Cl_lumen, Cl_stroma, Hstroma, pHstroma=y
    # pHlumen, Dy, pmf, Klumen, Kstroma, ATP_made, ATP_pool, ADP_pool, Cl_lumen, Cl_stroma,Hstroma, pHstroma = y

    #pmf计算
    # pmf=Dy + 0.06*(pHstroma-pHlumen) 

    #KEA3
    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    v_KEA = k_KEA*(Hlumen*Kstroma -  Hstroma*Klumen)
   
    #V_K 可能要加调控
    K_deltaG=-0.06*np.log10(Kstroma/Klumen) + Dy
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2
    
    #VCCN1
    driving_force_Cl = 0.06* np.log10(Cl_stroma/Cl_lumen) + Dy
    v_VCCN1 = k_VCCN1 * computer.Cl_flux_relative(driving_force_Cl) * (Cl_stroma + Cl_lumen)/2

    #CLCE
    v_CLCE =  k_CLCE*(driving_force_Cl*2+pmf)*(Cl_stroma + Cl_lumen)*(Hlumen+Hstroma)/4   

    #dK+/dt
    net_Klumen =  v_KEA - v_K_channel        
    dKlumen = net_Klumen*lumen_protons_per_turnover  
    dKstroma=0

    #dCl-/dt
    net_Cl_lumen_in = v_VCCN1 + 2*v_CLCE
    dCl_lumen = net_Cl_lumen_in * lumen_protons_per_turnover
    dCl_stroma = -0.1*dCl_lumen

    # dATP
    # activity = computer.ATP_synthase_actvt(t, T_ATP)
    activity = computer.ATP_synthase_actvt(t, 200)
    d_protons_to_ATP = computer.Vproton_pmf_actvt(pmf, activity, ATP_synthase_max_turnover, n)
    # d_ATP_made=d_protons_to_ATP/n 
    # d_ATP_consumed = d_ATP_made
    # dATP_pool= d_ATP_made - d_ATP_consumed
    # dADP_pool= - dATP_pool

    #dpHlumen
    # d_H_ATP_or_passive = computer.V_H_light(light_per_L, d_protons_to_ATP, pmf, Hlumen)      
    d_H_ATP_or_passive = computer.V_H_light(d_protons_to_ATP, pmf, Hlumen)                         
    net_protons_in =  - d_H_ATP_or_passive
    dHin = (net_protons_in - v_KEA - v_CLCE)*lumen_protons_per_turnover
    dpHlumen= -1*dHin / buffering_capacity 

    #dpHstroma
    dHstroma = 0
    dpHstroma = -1*dHstroma / buffering_capacity

    #dDy
    delta_charges= - v_K_channel - v_VCCN1-3*v_CLCE 
    dDy=delta_charges*Volts_per_charge

    #dpmf
    dpmf= 0.06* dpHlumen + dDy

    return [dpHlumen, dDy, dpmf, dKlumen, dKstroma, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]


# 模拟函数
def simulate(initial_states, t_span, params):
    sol = solve_ivp(f, t_span, initial_states, args=params, method='BDF', dense_output=True)
    return sol

# 数据处理与绘图
def plot_results(sol, labels):
    """绘制模拟结果"""
    time = sol.t
    results = sol.y

    plt.figure(figsize=(12, 8))
    for i, label in enumerate(labels):
        plt.plot(time, results[i], label=label)
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration/Gradient')
    plt.title('Simulation Results')
    plt.legend()
    plt.show()



# 设置初始条件和参数
initial_states = [7.8, 0.0, 0.0, 0.1, 0.1, 0.04, 0.04, 0.0, 7.8]
params = (0.000587, 200.0, 0.047, 150, 4.666, 0.03, 2500000, 12, 800000) 
t_span = (0, 1200)  # 模拟时间范围

# 运行模拟并绘制结果
species_labels = ['pHlumen', 'Dy', 'pmf', 'Klumen', 'Kstroma', 'Cl_lumen', 'Cl_stroma', 'Hstroma', 'pHstroma']
sol = simulate(initial_states, t_span, params)
plot_results(sol, species_labels)