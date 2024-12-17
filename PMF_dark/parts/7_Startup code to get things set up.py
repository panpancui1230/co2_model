# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:54:11 2024

@author: LENOVO
"""

   
#*******************************************************************************
#*******************************************************************************
#                    Startup code to get things set up. 
#*******************************************************************************
#********************************************************************************
#global FREQUENCY

def do_stuff(LIGHT):
    """
    to simulate a genotype, define an empty dictionary, then run sim_a_gtype()
    """
    print(LIGHT)
    WT = {}
    sim_a_gtype(WT, 'WT', LIGHT)
    # clce2 = {}
    # sim_a_gtype(clce2, 'clce2', 500)
    kea3 ={}
    sim_a_gtype(kea3, 'kea3', LIGHT)
    # vccn1 = {}
    # sim_a_gtype(vccn1,'vccn1', LIGHT)
    # # cckk ={}
    # # sim_a_gtype(cckk, 'clce2kea3', 500)
    # # ccvv ={}
    # # sim_a_gtype(ccvv,'clce2vccn1', 500)
    # v1k3 = {}
    # sim_a_gtype(v1k3, 'vccn1kea3', LIGHT)
    # vck ={}
    # sim_a_gtype(vck, 'vccn1clce2kea3', 500)
    
    
    """
    to run a simulation at 500 uE, one needs to ajust the T for ATP_synthase_actvt(t)
    function, default is 165s for 100uE, change it to 60s for 500uE
    """
    
    #####this following code saves NPQ difs between mutants and WT, simulated####
    ##### it can be easily modified to save other difs#####
    # mutant_list = [ kea3, vccn1, v1k3]
    # mutant_strs = [ 'kea3','vccn1','v1k3']
    # df_list = []
    time_min = WT['time_axis']/60
    idx = np.argwhere(time_min == 2)[0][0]
    
    delta_NPQ = kea3['NPQ']-WT['NPQ']
    delta_LEF = kea3['LEF']-WT['LEF']
    
    # df_list.append(time_min)
    df_ = {}
    # for idx, a_mutant in enumerate(mutant_list):
    #     delta_NPQ = a_mutant['NPQ'] - WT['NPQ']
    #     #df_NPQ_list = df_NPQ.tolist()
    df_['kea3_dNPQ'] = delta_NPQ
    df_['kea3_dLEF'] = delta_LEF
    df_['WT_NPQ'] = WT['NPQ']
    df_['WT_LEF'] = WT['LEF']
    #     df_list.append(df_NPQ)
    fig = plt.figure(num=3, figsize=(5,4), dpi=200)
    plt.plot(time_min[1:],delta_NPQ[1:],label = '∆NPQ: kea3 - WT')
    plt.legend()
    plt.show()
    plt.close()
    fig = plt.figure(num=3, figsize=(5,4), dpi=200)
    plt.plot(time_min[1:],delta_LEF[1:],label = '∆LEF: kea3 - WT')
    plt.legend()
    plt.show()
    plt.close()
    pdindex =pd.Index(WT['time_axis'], name = 'time/s')
    df_NPQ = pd.DataFrame(df_, index = pdindex)
    # #list_to_csv(df_list,str(1/FREQUENCY)+'df_NPQ_sin_500_simulated.csv')
    df_NPQ.to_csv('delta_NPQ_LEF'+str(LIGHT)+'_uE_simulated.csv')
    plt.show()
    plt.close()
    
    #a =   df_NPQ[df_NPQ.index>= df_NPQ.index[-1]-1/FREQUENCY].max()\
    #    - df_NPQ[df_NPQ.index>= df_NPQ.index[-1]-1/FREQUENCY].min()        
    return (delta_NPQ[idx], delta_LEF[idx],\
            delta_NPQ[idx]/WT['NPQ'][idx], delta_LEF[idx]/WT['LEF'][idx])

global FREQUENCY, LIGHT, T_ATP
FREQUENCY = 1/60
result_dict = {}
light_T = [(50, 200), (100, 165), (250, 100), (500, 60), (1000, 40)]
for LIGHT, T_ATP in light_T:
    delta = do_stuff(LIGHT)
    result_dict[LIGHT] = delta
col_list = ['dNPQ_2min', 'dLEF_2min', 'dNPQ_rel', 'dLEF_rel']    
column = {}
for i, col in enumerate(col_list):
    column[i] = col
result_df = pd.DataFrame(result_dict).T
result_df.rename(columns = column, inplace= True)
#result_df.to_csv('/Users/LIMeng/Desktop/KEA3_Qs/2min_delta_NPQ_LEF_kea3_abs_relative.csv')
