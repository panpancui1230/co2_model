import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd


from painter import Plotting
from utils import standard_constants
from utils import standard_initial_states
from utils import sim_states
from utils import sim_constants

from calc import block

from sun_sim import sunshine
# labels for the results of odeint(f, ... )
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

# 初始状态
initial_states = [7.8, 0.0, 0.0, 0.1, 0.1, 0.04, 0.04, 0.0, 7.8]
# initial_states = [6.5, 50.0, 1.0, 0.5, 0.5, 0.01, 0.01, 1.0, 8.0]
# 参数
params = [0.000587, 200.0, 0.047, 150, 4.666, 0.03, 2500000, 12, 800000]

# 时间范围
t_end = 12000
# time_eval = np.linspace(0, t_end, t_end + 1)  # 每秒一个点
time_eval = np.linspace(0, t_end, t_end * 10 + 1)  # 12000点，每秒10个点

def f(t, y, lumen_protons_per_turnover, ATP_synthase_max_turnover, Volts_per_charge, perm_K,
      n, buffering_capacity, k_KEA, k_VCCN1, k_CLCE):
    computer = block()

    pHlumen, Dy, pmf, Klumen, Kstroma, Cl_lumen, Cl_stroma, Hstroma, pHstroma=y

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

# 求解 ODE
solution = solve_ivp(f, [0, t_end], initial_states, t_eval=time_eval, args=tuple(params))

# 提取结果
time = solution.t
results = solution.y

# 准备数据以导出为 CSV
output_data = {
    'Time (s)': time,
    'pH lumen': results[0],
    'Dy': results[1],
    'pmf': results[2],
    'Klumen': results[3],
    'Kstroma': results[4],
    'Cl lumen': results[5],
    'Cl stroma': results[6],
    'Hstroma': results[7],
    'pH stroma': results[8]
}

df = pd.DataFrame(output_data)

# 保存为 CSV 文件
output_file = "./simulation_results.csv"
df.to_csv(output_file, index=False)

print(f"Simulation data saved to {output_file}")