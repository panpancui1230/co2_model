import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import solve_ivp
from scipy.integrate import cumulative_trapezoid
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

warnings.filterwarnings("ignore")
max_light_change=1
points_per_segment=1000


#Function f calculates the changes in state for the entire systems
def f(t, y, ratio_absorb,PSII_content_per_leaf,pKreg, max_PSII, kQA, max_b6f, lumen_protons_per_turnover, PAR, ATP_synthase_max_turnover, 
    PSII_antenna_size, Volts_per_charge, perm_K, n, Em7_PQH2, Em7_PC,Em_Fd, PSI_antenna_size, 
    buffering_capacity, VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pKPsbS, max_NPQ, k_recomb, k_PC_to_P700, 
    triplet_yield, triplet_to_singletO2_yield, fraction_pH_effect, k_Fd_to_NADP, k_CBC, k_KEA, k_VCCN1, k_CLCE, k_NDH): 
    sun = sunshine()
    #The following are holders for paramters for testing internal functions of f
    PAR = sun.light(t, 1200, LIGHT, FREQUENCY, 900, 100)
    light_per_L=ratio_absorb * PAR/PSII_content_per_leaf


    computer = block()
    QA, QAm, PQ, PQH2, Hin, pHlumen, Dy, pmf, deltaGatp, Klumen, Kstroma, ATP_made,\
    PC_ox, PC_red, P700_ox, P700_red, Z, V, NPQ, singletO2, Phi2, LEF, Fd_ox, Fd_red,\
    ATP_pool, ADP_pool, NADPH_pool, NADP_pool,Cl_lumen, Cl_stroma, Hstroma, pHstroma =y
    
    PSII_recombination_v=computer.recombination_with_pH_effects(k_recomb, QAm, Dy, pHlumen, fraction_pH_effect)
    dsingletO2=PSII_recombination_v*triplet_yield*triplet_to_singletO2_yield

    #calculate pmf from Dy and deltapH 
    pmf=Dy + 0.06*(pHstroma-pHlumen)

    # Phi2=computer.Calc_Phi2(QA, NPQ) #I use the current' value of NPQ. I then calculate the difference below 
    Phi2=computer.Calc_Phi2(QA,QAm,NPQ)
    PSII_charge_separations=PSII_antenna_size*light_per_L * Phi2
    
    
    Keq_QA_PQ=200
    
    #calculate the changes in QA redox state based on the number of charge separations and equilibration with 
    #the PQ pool
    dQAm = PSII_charge_separations  + PQH2*QA*kQA/Keq_QA_PQ  - QAm * PQ * kQA - PSII_recombination_v
    dQA = -1*dQAm

    b6f_content=0.433 #Journal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014

    v_b6f=computer.calc_v_b6f(max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf)
    
    v_NDH = computer.calc_v_NDH(Em_Fd, Em7_PQH2, pHstroma, pmf, k_NDH, Fd_red, Fd_ox, PQ, PQH2)
    d_Hlumen_NDH = v_NDH*2 #change in lumen protons
    d_charge_NDH = d_Hlumen_NDH # change in charges across the membrane
    #d_Hstroma_NDH = v_NDH*3 # change in stroma protons
    
    ##PGR regulation, attempted
    PGR_vmax = 0#It seems this function does not impact the kinetics much.
    v_PGR = computer.calc_v_PGR(PGR_vmax, Fd_red, PQ, PQH2)

    #calculate the change in PQH2 redox state considering the following:
    #PQ + QAm --> PQH2 + QA ; PQH2 + b6f --> PQ    
    PSI_charge_separations= P700_red * light_per_L * PSI_antenna_size * Fd_ox

    dPQH2 = (QAm * PQ * kQA + v_NDH + v_PGR - v_b6f - PQH2*QA*kQA/Keq_QA_PQ)*0.5 
    dPQ = -1*dPQH2
    
    
    #P700 reactions
    d_P700_ox = PSI_charge_separations - PC_red * k_PC_to_P700 * P700_ox
    d_P700_red=-1*d_P700_ox
    
    #PC reactions:
    d_PC_ox = PC_red * k_PC_to_P700 * P700_ox - v_b6f
    d_PC_red = -1*d_PC_ox
    
    #Mehler reaction, V_me = kme * [O2]*Fd_red/(Fd_red+Fd_ox), Hui Lyu and Dusan Lazar modeling...
    V_me = 4*0.000265*Fd_red/(Fd_red+Fd_ox)
    dFd_red=PSI_charge_separations - k_Fd_to_NADP*Fd_red*NADP_pool - v_NDH - v_PGR -V_me
    dFd_ox=-1*dFd_red

    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    activity = computer.ATP_synthase_actvt(t, T_ATP)
    
    d_protons_to_ATP = computer.Vproton_pmf_actvt(pmf, activity, ATP_synthase_max_turnover, n)
    d_H_ATP_or_passive = computer.V_H_light(light_per_L, d_protons_to_ATP, pmf, Hlumen)                              

        
    d_ATP_made=d_protons_to_ATP/n                                        

    NADPH_CBC = k_CBC*(1.0-np.exp(-t/900))*(np.log(NADPH_pool/NADP_pool)-np.log(1.25))/(np.log(3.5/1.25))#calc_CBC_NADPH(k_CBC, t, d_ATP_made)
    #this number in "np.exp(-t/600)" is important, which impacts the shape of the curves
    dNADPH_pool=0.5 * k_Fd_to_NADP*NADP_pool*Fd_red - NADPH_CBC
    dNADP_pool=-1*dNADPH_pool
    
    dLEF=k_Fd_to_NADP*NADP_pool*Fd_red
    
    d_ATP_consumed = d_ATP_made#NADPH_CBC*5/3 + (ATP_pool/(ADP_pool+ATP_pool)-0.5)*1.2#ATP_pool*(ATP_pool/ADP_pool-1)  
    d_protons_from_PSII = PSII_charge_separations - PSII_recombination_v

    #calculate the contributions to Dy from PSII
    charges_from_PSII = PSII_charge_separations - PSII_recombination_v
    

    d_protons_from_b6f = v_b6f*2 #two protons per electron transferred from PQH2 to PC


    charges_from_b6f = v_b6f
     
    #add up the changes in protons delivered to lumen
    #note: net_protons_in is the total number of protons input into the lumen, including both free and bound.
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

    dHstroma = 0#(net_protons_stroma + v_KEA + v_CLCE)*lumen_protons_per_turnover/10
    #Assuming the volume of stroma is ten times as that of lumen
    dpHstroma = -1*dHstroma / buffering_capacity
    delta_charges=charges_from_PSII+PSI_charge_separations + charges_from_b6f \
                    + d_charge_NDH - v_K_channel - d_H_ATP_or_passive - v_VCCN1-3*v_CLCE
    
    dDy=delta_charges*Volts_per_charge
    dpmf= 0.06* dpHlumen + dDy


    dATP_pool= d_ATP_made - d_ATP_consumed
    dADP_pool= - dATP_pool
    
    dZ, dV = computer.calc_v_VDE(VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pHlumen, V, Z)

    
    new_PsbS_H = computer.calc_PsbS_Protonation(pKPsbS, pHlumen + dpHlumen)
    new_Z=Z+dZ
    
    new_NPQ=0.4*max_NPQ*new_PsbS_H*new_Z+0.5*max_NPQ*new_PsbS_H+0.1*max_NPQ*new_Z
    dNPQ=new_NPQ-NPQ #new_PsbS_H-PsbS_H
    dPhi2=0 #

    ddeltaGatp = 0
    return [dQA, dQAm, dPQ, dPQH2, dHin, dpHlumen, dDy, dpmf, ddeltaGatp, dKlumen, dKstroma, 
            d_ATP_made, d_PC_ox, d_PC_red, d_P700_ox, d_P700_red, dZ, dV, dNPQ, dsingletO2, dPhi2, dLEF, 
            dFd_ox, dFd_red,  dATP_pool, dADP_pool, dNADPH_pool,dNADP_pool, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]
            
def sim(K, initial_states, pulse_times_and_light, max_light_change=1, points_per_segment=1000, **keyword_parameters):    
    
    sun = sunshine()
    if ('dark_equilibration' in keyword_parameters):
        equibrate_time= keyword_parameters['dark_equilibration']
        use_initial_states=dark_equibration(initial_states, K, equibrate_time)
    else:
        use_initial_states=initial_states
    sub_arrays_time_and_light= sun.optimized_time_split(pulse_times_and_light, 
                                max_light_change, points_per_segment) 

    #first make the constants_set and trace_times for the split segments

    constants_set_and_trace_times=sun.make_variable_light_constants_set_and_trace_times(K, sub_arrays_time_and_light)

    #next, do the simulation
    output=do_complete_sim(use_initial_states, constants_set_and_trace_times, K)
    return(output, use_initial_states)

def sim_ivp(K, initial_states, t_end):    
    output=do_complete_sim(initial_states, t_end, K)
    return output

def do_complete_sim(y00, t_end, Kx):

    computer = block()
    sun = sunshine()
    
    y00[11]=0 #set ATP_made to zero
    
    # prepare a dictionary with empty arrays to store output
    output={}
    for label in species_labels:
        output[label] = []
    
        
    soln = solve_ivp(f, [0, t_end], y00, args=Kx.as_tuple(), method = 'BDF', \
                     t_eval = np.linspace(0, t_end, 10*t_end+1), max_step = 5)


    time_axis = soln.t
            
        #append a set of computed constants to the output arrays
    for index, label in enumerate( species_labels ):
        output[label] = np.append( output[label], soln.y[index,:] )

    # save the results in case we want to start another simulation that starts where this one left off
    end_state = list(soln.y[:,-1])

    Dy = output['Dy']
    pHlumen = output['pHlumen']
    pHstroma = output['pHstroma']
    pmf_total= Dy + ((pHstroma-pHlumen)*.06)

    Phi2_array=[] #contains the calculated Phi2 results
    QA = output['QA']
    QAm=output['QAm']
    NPQ_array = output['NPQ_array']

    for i in range(len(QA)):
        # Phi2_array.append(computer.Calc_Phi2(QA[i], NPQ_array[i]))
        Phi2_array.append(computer.Calc_Phi2(QA[i],QAm[i],NPQ_array[i]))
        
    #calculate tPhiNO and PhiNPQ
    # using the Calc_PhiNO_PhiNPQ function.
    PhiNPQ_array=[]
    PhiNO_array=[]
    for i in range(len(QA)):
        PhiNO, PhiNPQ=computer.Calc_PhiNO_PhiNPQ(Phi2_array[i], QA[i], NPQ_array[i])
        PhiNPQ_array.append(PhiNPQ)
        PhiNO_array.append(PhiNO)
    output['PhiNPQ']=PhiNPQ_array
    output['PhiNO']=PhiNO_array

    #Set up an array to contain the light curve (the PAR values), 
    light_curve=[]    
    for a_t in time_axis:
        light_curve.append(sun.light(a_t, 1200, LIGHT, FREQUENCY, 900, 100))

    # compute LEF array from Phi2 and light (does not consider recombination!)
    LEF_array_from_Phi2=[]
    PSII_cross_section=0.425
    for i in range(0,len(Phi2_array)):
        LEF_array_from_Phi2.append(light_curve[i]*Phi2_array[i]*PSII_cross_section)
    ###how about PSII_cross_section = 0.5 and then times leaf absorbance of 0.85?


    # calculate singletO2_rate
    singletO2_array = output['singletO2_array']
    singletO2_rate=[]
    singletO2_rate.append(0)
    for i in range(1,len(singletO2_array)):
        so2r=(singletO2_array[i]-singletO2_array[i-1])/(time_axis[i]-time_axis[i-1])
        singletO2_rate.append(so2r)
        
    # in singletO2_rate, get rid of nans when the delta_t was zero; replace with previous value
    for i in range(1,len(singletO2_rate)):
        if np.isnan(singletO2_rate[i]):
            singletO2_rate[i]=singletO2_rate[i-1]

    # compute delta_pH and delta_pH_V
    delta_pH=[]
    delta_pH_V=[]
    #fraction_Dy=[]
    for i in range(0,len(pHlumen)):
        dpH=pHstroma[i]-pHlumen[i]
        dpH_V=dpH*.06
        delta_pH.append(dpH)
        delta_pH_V.append(dpH_V)

    # before returning output, append the computed data
    # the output already includes all the results directly from odeint(f, ... )
    
    output['delta_pH']=delta_pH
    output['delta_pH_V']=delta_pH_V   
     
    output['pmf']=pmf_total

    output['delta_pH_offset']=delta_pH-delta_pH[0]
    output['delta_pH_V_offset']=delta_pH_V-delta_pH_V[0]    
    output['pmf_offset']=pmf_total-pmf_total[0]
    output['Dy_offset']=Dy-Dy[0]
    
    #output['deltaGatp']=deltaGatp
    output['pmf_total']=pmf_total
    output['singletO2_rate']=singletO2_rate
    output['time_axis']=time_axis
    output['time_axis_min']=time_axis/60
    output['time_axis_h']=time_axis/3600
    output['end_state'] = end_state
    output['light_curve'] = light_curve
    integrated_light=[]
    il_cum=0.0
    integrated_light.append(il_cum)
    for indexl in range(1, len(light_curve)):
        il_cum=il_cum+light_curve[indexl]*(time_axis[indexl]-time_axis[indexl-1])
        integrated_light.append(il_cum)
    output['integrated_light']=integrated_light
    output['fraction_Dy']=output['Dy']/output['pmf']
    output['fraction_DpH']=1-output['fraction_Dy']

    
    output['Z']=output['Z_array'] 
    output['V']=output['V_array'] 
    output['NPQ']=NPQ_array 
    output['singletO2']=singletO2_array 
    output['Phi2'] = Phi2_array
    
    output['LEF'] = np.array(LEF_array_from_Phi2) #output['LEF_array']    
    output['LEF_productive']=[] #output['LEF_array']
    output['LEF_productive'].append(0)

    for i in range(1,len(output['LEF_array'])):
        temp=(output['LEF_array'][i]-output['LEF_array'][i-1])/(time_axis[i]-time_axis[i-1])
        output['LEF_productive'].append(temp)

    #calculate the electron flow to NADPH
    output['LEF_to_NADPH']=[0]
    #output['LEF_to_NADPH'].append(0)
    for i in range(1,len(output['LEF_array'])):
        temp=(output['LEF_array'][i]-output['LEF_array'][i-1])/(time_axis[i]-time_axis[i-1])
        output['LEF_to_NADPH'].append(temp)
    output['LEF_to_NADPH']=np.array(output['LEF_to_NADPH']) #convert to np array so we can do calculations below
    
    LEF_cumulative=[]
    LEF_cumulative.append(0)
    LEF_cum=0
    
    ATP_rate=[0]
    for i in range(1, len(output['LEF'])):
        LEF_cum=LEF_cum+(output['LEF'][i] * (output['time_axis'][i]-output['time_axis'][i-1]))
        LEF_cumulative.append(LEF_cum)
        q1=np.array(output['ATP_made'][i])-np.array(output['ATP_made'][i-1])
        q2=np.array(output['time_axis'][i])-np.array(output['time_axis'][i-1])
        delta_ATP=(q1/q2)
        ATP_rate.append(delta_ATP)
    normalized_LEF_cumulative=LEF_cumulative/LEF_cumulative[-1]
    output['LEF_cumulative']=output['LEF_array'] #LEF_cumulative
    output['normalized_LEF_cumulative']=normalized_LEF_cumulative
    output['ATP_rate']=np.array(ATP_rate)
    NADPH=np.array(output['LEF'], dtype=float)/2
    output['NADPH']=NADPH
   
    output['PsbS_protonated']=[]
    output['b6f_control']=[]
    for pH in output['pHlumen']:
        output['PsbS_protonated'].append(computer.calc_PsbS_Protonation(Kx.pKPsbS, pH))
        output['b6f_control'].append(computer.calc_k_b6f(Kx.max_b6f,1, pH, Kx.pKreg))
        
    Fd_rate=[0]
    for index in range(1,len(output['Fd_red'])):
        Fd_rate.append((output['Fd_red'][index]-output['Fd_red'][index-1])/(output['time_axis'][index]-output['time_axis'][index-1]))
    output['Fd_rate']=np.array(Fd_rate)
    output['ATP/NADPH']= 2*output['ATP_rate']/(output['Fd_rate'])
    
    K_flux=[0] #start with zero because we will be calculating the derivative
    for i in range(1,len(output['Klumen'])):
        K_flux.append((output['Klumen'][i-1]-output['Klumen'][i])/(output['time_axis'][i]-output['time_axis'][i-1]))

    output['K_flux']=np.array(K_flux)
    output['K_flux']=output['K_flux']/Kx.lumen_protons_per_turnover
    
    for i in range(len(output['ATP_rate'])):
        if np.isnan(output['LEF_to_NADPH'][i]):
            output['LEF_to_NADPH'][i]=output['LEF_to_NADPH'][i-1]
        if np.isnan(output['ATP_rate'][i]):
            output['ATP_rate'][i]=output['ATP_rate'][i-1]
        if np.isnan(output['K_flux'][i]):
            output['K_flux'][i]=output['K_flux'][i-1]

    # calculate the deficit in ATP/NADPH and store it in output['deficit']
    output['deficit']=(output['LEF_to_NADPH']*(3.0/Kx.n)-output['ATP_rate'])  
    output['deficit_int']=integrate.cumulative_trapezoid(output['deficit'], output['time_axis'], initial=0)
    output['fract_deficit']=output['deficit_int']/output['LEF_to_NADPH']

    return(output)

def dark_equibration(y_initial, Kx, total_duration, **keyword_parameters): 
    #make a sin wave with zero amplitude

    sun = sunshine()
    light_frequency=1/total_duration
    points_per_second=10
    max_PAR=0
    dark_time_light_profile=sun.generate_sin_wave(total_duration, max_PAR, light_frequency, points_per_second)
    max_light_change=10
    points_per_segment=100
    optimized_dark_sub_arrays= sun.optimized_time_split(dark_time_light_profile, 
        max_light_change, points_per_segment) 

    constants_set_and_times = sun.make_variable_light_constants_set_and_trace_times(Kx, optimized_dark_sub_arrays)

    output=do_complete_sim(y_initial, constants_set_and_times, Kx)
    dark_equilibrated_initial_y=output['end_state']

    if ('return_kinetics' in keyword_parameters) and keyword_parameters['return_kinetics']==True:
        return(dark_equilibrated_initial_y, output)
    else:
        return(dark_equilibrated_initial_y)

    
# global_painter = Plotting()
# plot_results={}
# plot_results['pmf_params']=global_painter.plot_pmf_params
# plot_results['pmf_params_offset']=global_painter.plot_pmf_params_offset
# plot_results['K_and_parsing']=global_painter.plot_K_and_parsing
# plot_results['plot_QAm_and_singletO2']=global_painter.plot_QAm_and_singletO2
# plot_results['plot_cum_LEF_singetO2']=global_painter.plot_cum_LEF_singetO2
# plot_results['b6f_and_balance'] = global_painter.b6f_and_balance
    

class ListTable(list):
    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")
            for col in row:
                html.append("<td>{0}</td>".format(col))
            
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)
        
#display only the constants that are different
# def Changed_Constants_Table(table_title, original_values, Kxx):
#     table = ListTable()
#     table.append(['Changed Parameter', 'Old Value', 'New Value']) #, 'Short Description'])
#     Kdict=Kxx.as_dictionary()
#     original_values_dict=original_values.as_dictionary()

#     for key in list(Kdict.keys()):
#         if Kdict[key] == original_values_dict[key]:
#             pass
#         else:
#             table.append([key, original_values_dict[key], Kdict[key]]) #, Ksumm[key]])
#     print(table_title)
#     display(table)

#PROCESS A GENOTYPE AND SAVE ITS SIMULATED VALUES
def process_a_gtype(gtype_dict, parameter_list, out_dict, gtype='a_genotype'):
    gtype_df = pd.DataFrame([])
    for para in parameter_list:
        gtype_dict[para] = out_dict[para]#store in dictionary for further calculation
        gtype_df[para] = [np.mean(out_dict[para][-20:])]
    
    return gtype_df


def sim_a_gtype(gtype_dict, idx=0, light = 100):  
    parameters_of_interest = ['time_axis','NPQ','Phi2','LEF','qL','Z','V',\
                          'pmf','Dy','pHlumen','fraction_Dy','fraction_DpH',\
                          'Klumen','Cl_lumen','Cl_stroma']

    initial_sim_states=sim_states()
    # initial_sim_state_list=initial_sim_states.as_list()
  
    Kx_initial=sim_constants() 
       

    constants_dict={}
    k_CBC_light = 60 * (light/(light+250))#this needs change with different light intensity    

    output_dict={}
    on = str(idx)
    Kx=sim_constants()
    # if 'clce2' in gtype:
    #     Kx.k_CLCE = 0
    # if 'kea3' in gtype:
    #     Kx.k_KEA =0
    # if 'vccn1' in gtype:
    #     Kx.k_VCCN1 =0
    Kx.k_CBC = k_CBC_light
    csv_file='./data/constants.csv'
    data_df = pd.read_csv(csv_file) # get the changing constants dataframe
    varying = data_df.loc[idx, data_df.columns[:3]]  # take row index, and first 3 columns
    
    

    # now change the constants
    Kx.ratio_absorb = varying[0]
    Kx.PSII_content_per_leaf = varying[1]
    Kx.PSII_antenna_size = varying[2]

    new_scv_file = './data/initial_states.csv'
    data_is = pd.read_csv(new_scv_file)
    # new_verifying = data_is.loc[idx, data_is.columns[:1]]
    # Kx.P700_red_initial = new_verifying[0]
    # initial_sim_state_list[0]=new_verifying[0]

    # initial_sim_state_list[0] = data_is.loc[idx, 'P700_red_initial']

    initial_sim_state_list['P700_red_initial'] = data_is.loc[idx, 'P700_red_initial']
    initial_sim_state_list['QA_content_initial'] = data_is.loc[idx, 'QA_content_initial']
    initial_sim_state_list=initial_sim_states.as_list()
    # initial_sim_state_list[15] = data_is.loc[idx, 'P700_red_initial']
    # print(initial_sim_state_list[15])
    
    # Kx.PSII_content_per_leaf = new_verifying[1]
    # Kx.PSII_antenna_size = new_verifying[2]


    
    # 输出当前的 ratio_absorb、PSII_content_per_leaf 和 PSII_antenna_size 值
    # print(f"Current parameters for idx {idx}:")
    # print(f"ratio_absorb: {Kx.ratio_absorb}")
    # print(f"PSII_content_per_leaf: {Kx.PSII_content_per_leaf}")
    # print(f"PSII_antenna_size: {Kx.PSII_antenna_size}")
    # print(f"P700_red_initial:{initial_sim_state_list[0]}")
    # print(f"P700_red_initial")

    constants_dict[on]=Kx #store constants in constants_dict

    output_dict=sim_ivp(Kx, initial_sim_state_list, 1200)
    # Changed_Constants_Table('Change Constants', Kx_initial, Kx)
    # output_dict['qL'] = 1-output_dict['QAm']

    output_dict['qL'] = output_dict['QA']/(output_dict['QA']+output_dict['QAm'])


    # paint = Plotting()
    # paint.plot_interesting_stuff(str(idx), output_dict)

    # plot_interesting_stuff(gtype, output_dict)
    # process_a_gtype(gtype_dict,parameters_of_interest, output_dict, str(idx)+'_'+str(light)+'uE')   
    return process_a_gtype(gtype_dict,parameters_of_interest,output_dict,str(idx)+'_'+str(light)+'uE')  



# 修改 do_stuff 函数，累积 1000 个仿真结果并保存到单个 CSV 文件
# def do_stuff(LIGHT):
#     print(f"Running simulations for LIGHT intensity: {LIGHT}")
#     combined_df = pd.DataFrame()  # 用于存储所有 idx 的 gtype_df 数据
#     for idx in range(840):  # 循环从索引 0 到 999
#         gtype_dict = {}
#         gtype_df = sim_a_gtype(gtype_dict, idx=idx, light=LIGHT)  # 获取单个仿真结果的 gtype_df
#         gtype_df['idx'] = idx  # 添加 idx 列以标识每个组合
#         combined_df = pd.concat([combined_df, gtype_df], ignore_index=True)  # 累积到 combined_df
    
#     # 保存所有 idx 的结果到一个 CSV 文件
#     combined_df.to_csv(f'./logs_cpp/combined_{LIGHT}_simulated_cpp_6.0.csv', index=False)
#     print(f"All results for LIGHT {LIGHT} saved to combined_{LIGHT}_simulated_cpp_6.0.csv")


# 修改 do_stuff 函数，累积 1000 个仿真结果并保存到单个 CSV 文件
def do_stuff(LIGHT):
    print(f"Running simulations for LIGHT intensity: {LIGHT}")
    combined_df = pd.DataFrame()  # 用于存储所有 idx 的 gtype_df 数据
    # error_indices = []
    for idx in range(840):  # 循环从索引 0 到 839
        try:
            gtype_dict = {}
            gtype_df = sim_a_gtype(gtype_dict, idx=idx, light=LIGHT)  # 获取单个仿真结果的 gtype_df
            gtype_df['idx'] = idx  # 添加 idx 列以标识每个组合
            combined_df = pd.concat([combined_df, gtype_df], ignore_index=True)  # 累积到 combined_df
        except ValueError as e:
            print(f"ValueError at idx {idx}: {e}")
            # error_indices.append(idx)
        # except Exception as e:
        #     print(f"An error occurred at idx {idx}: {e}")
        #     error_indices.append(idx)

    # 保存所有 idx 的结果到一个 CSV 文件
    combined_df.to_csv(f'./logs_QA/combined_{LIGHT}_simulated_1.csv', index=False)
    print(f"All results for LIGHT {LIGHT} saved to combined_{LIGHT}_simulated_1.csv")

global FREQUENCY, LIGHT, T_ATP
FREQUENCY = 1/60
result_dict = {}
light_T = [(100, 165)]
# light_T = [(500, 60)]
# light_T = [(100, 165),(500, 60)]
# light_T = [(50, 200), (100, 165), (250, 100), (500, 60), (1000, 40)]
for LIGHT, T_ATP in light_T:
    do_stuff(LIGHT)
