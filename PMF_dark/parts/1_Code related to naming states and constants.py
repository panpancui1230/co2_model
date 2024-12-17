# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:27:40 2024

@author: LENOVO
"""


# -*- coding: utf-8 -*-

#*******************************************************************************
#*******************************************************************************
#                       Importing required libraries                           *
#*******************************************************************************
#*******************************************************************************


# -*- coding: utf-8 -*-
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


#*******************************************************************************
#*******************************************************************************
#                   Code related to naming states and constants                *
#*******************************************************************************
#*******************************************************************************

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


#Location for the PDF files describing the parameters
#This can be over-written during the actual run
#PDF_file_location='/Users/davidkramer/Dropbox/Data/DF_ECS/Params.png/' 

#In several places the code sometimes returns Nans resulting from divisions by 
#zero. This code supresses the warnings to de-clutter the display.

import warnings
warnings.filterwarnings("ignore")




#*******************************************************************************
#*******************************************************************************
#                   Code related to generating light curves                    *
#*******************************************************************************
#*******************************************************************************

#the following two parameters are used to generate the time and light curves
#for the simulations. It is needed to prevent very abrupt changes in conditons
#that would make the simulaitons too stiff for odeint. The default values are 
#typically OK, but can be adjusted if odeint produces garbage.

max_light_change=1
points_per_segment=1000


#the following code sets up a light regime based on sin waves. The wave can 
#be either a sin or a square wave.

#total_duration=the duraction in time units of the complete wave form
#time_units are the time units, either  'seconds', 'minutes' or 'hours'.
#max_PAR is the maximum light intensity
#wave_form indicates if the wave is either a 'sin' or a 'square' wave
#light_frequency is the frequencgy in Hz for the waveform
# e.g. light_frequency=1/(60*60) will result in a one-hour duration wave




#number_segments is the number of segments to split the sequence. The time_slice will be adjusted for each segment
#to keep the changes under a certain value.   

def generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity, pulse_duration, pulse_intensity, 
                                        recovery_duration, recovery_intensity, rise_time, time_units, point_frequency,
                                        repeat_cycles):
    pulse_times=[]
    pulse_light=[]

    if time_units=='seconds':
        time_div=1
    elif time_units=='minutes':
        time_div=60
    elif time_units=='hours':
        time_div=60*60

    baseline_duration=baseline_duration*time_div
    baseline_points=baseline_duration*point_frequency
    baseline_time=np.linspace(0, baseline_duration, baseline_points)
    baseline_intensity_array= np.linspace(baseline_intensity, baseline_intensity,baseline_points)
    pulse_times=np.append(pulse_times, baseline_time)
    pulse_light=np.append(pulse_light, baseline_intensity_array)

    riser_duration=rise_time*time_div
    riser_points=riser_duration*point_frequency
    riser_start_time = (baseline_points+1) / point_frequency
    riser_end_time = riser_start_time + riser_duration
    riser_time=np.linspace(riser_start_time, riser_end_time, riser_points)
    riser_light=np.linspace(baseline_intensity, pulse_intensity, riser_points)
    
    pulse_times=np.append(pulse_times, riser_time)
    pulse_light=np.append(pulse_light, riser_light)
    pulse_duration=pulse_duration*time_div
    pulse_points=pulse_duration*point_frequency
    pulse_start_time = (baseline_points + riser_points +1)/point_frequency
    pulse_end_time = pulse_start_time + pulse_duration
    pulse_time=np.linspace(pulse_start_time, pulse_end_time, pulse_points)
    pulse_light_array=np.linspace(pulse_intensity, pulse_intensity, pulse_points)
    pulse_times=np.append(pulse_times, pulse_time)
    pulse_light=np.append(pulse_light, pulse_light_array)
    
    falling_duration=rise_time*time_div
    falling_points=riser_duration*point_frequency
    falling_start_time = (baseline_points + riser_points + pulse_points + 1) / point_frequency
    falling_end_time = falling_start_time + falling_duration
    falling_time=np.linspace(falling_start_time, falling_end_time, falling_points)
    falling_light=np.linspace(pulse_intensity, recovery_intensity, falling_points)
    
    pulse_times=np.append(pulse_times, falling_time)
    pulse_light=np.append(pulse_light, falling_light)
    
    recovery_duration=recovery_duration*time_div
    recovery_points=recovery_duration*point_frequency
    recovery_start_time = (baseline_points + riser_points + pulse_points + falling_points + 1) / point_frequency
    recovery_end_time = recovery_start_time + recovery_duration
    recovery_time=np.linspace(recovery_start_time, recovery_end_time, recovery_points)
    recovery_light=np.linspace(recovery_intensity, recovery_intensity, recovery_points)

    pulse_times=np.append(pulse_times, recovery_time)
    pulse_light=np.append(pulse_light, recovery_light)
    pulse_times_seq=[]
    pulse_light_seq=[]
    
    for index in range(0,repeat_cycles):
        pulse_times_seq=np.append(pulse_times_seq, pulse_times + index * pulse_times[-1])
        pulse_light_seq=np.append(pulse_light_seq, pulse_light)
    return([pulse_times_seq, pulse_light_seq])



#smooths a trace using the simple boxcar algorithm. 'box_pts' is the box size and y is the trace. It assumes equally spaced data\n",
def smooth(yvals, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(yvals, box, mode='full')
    return y_smooth
    

    
# def make_waves():
#     # Generate a library of light patterns to use in the simulations.
#     light_pattern={}
    
#     #a single turnover flash  
    
#     baseline_duration=1 # 10 ms dark timein seconds
#     baseline_intensity=0 #dark baseline
#     pulse_duration=0.001 # 10 ms pulse of bright light 
#     pulse_intensity=3000 #pulse is 1000 units
#     recovery_duration = 10 #10 s recovery in dark
#     recovery_intensity=0 #recovery is dark
#     rise_time=.001 #100 ms for the light to rise
#     time_units='seconds' 
#     point_frequency=1000 #start with a frequency of 1000 points per subtrace
#     repeat_cycles=1 #do this once
#     wave=generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity,
#                         pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
#                         rise_time, time_units, point_frequency, repeat_cycles)
#     light_pattern['single_turnover_flash']=wave
    
    
#     #a single, 5-min square wave with peak intensity of 300 uE
#     baseline_duration=0 #in seconds
#     baseline_intensity=0 #dark baseline
#     pulse_duration=1200 #300 seconds pulse
#     pulse_intensity=100 #pulse is 1000 units
#     recovery_duration = 300 #100 seconds recovery
#     recovery_intensity=0 #recovery is dark
#     rise_time=1 #10 ms for the light to rise
#     time_units='seconds' 
#     point_frequency=100 #start with a frequency of 1000 points per subtrace
#     repeat_cycles=1 #do this once
#     wave=generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity,
#                         pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
#                         rise_time, time_units, point_frequency, repeat_cycles)
    
#     light_pattern['single_square_20_min_100_max']=wave

#     #a single, 20-min square wave with peak intensity of 225 uE/m2/s
#     baseline_duration=0 #in seconds
#     baseline_intensity=0 #dark baseline
#     pulse_duration=1200 #20 min light
#     pulse_intensity=225 #light intesntisy is 225 uE
#     recovery_duration = 300 #5 min recovery
#     recovery_intensity=0 #recovery is dark
#     rise_time=1 #10 ms for the light to rise
#     time_units='seconds' 
#     point_frequency=100 #start with a frequency of 1000 points per subtrace
#     repeat_cycles=1 #do this once
#     wave=generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity,
#                         pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
#                         rise_time, time_units, point_frequency, repeat_cycles)
#     light_pattern['single_square_20_min_225_max']=wave
    
#     #a single, 20-min square wave with peak intensity of 500 uE/m2/s
#     baseline_duration=0 #in seconds
#     baseline_intensity=0 #dark baseline
#     pulse_duration=1200 #20 min light
#     pulse_intensity=500 #light intesntisy is 500 uE
#     recovery_duration = 300 #5 min recovery
#     recovery_intensity=0 #recovery is dark
#     rise_time=1 #10 ms for the light to rise
#     time_units='seconds' 
#     point_frequency=100 #start with a frequency of 1000 points per subtrace
#     repeat_cycles=1 #do this once
#     wave=generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity,
#                         pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
#                         rise_time, time_units, point_frequency, repeat_cycles)
#     light_pattern['single_square_20_min_500_max']=wave
    
    
#     #a single, 20-min square wave with peak intensity of 500 uE/m2/s
#     baseline_duration=0 #in seconds
#     baseline_intensity=0 #dark baseline
#     pulse_duration=1200 #20 min light
#     pulse_intensity=0 #light intesntisy is 500 uE
#     recovery_duration = 300 #5 min recovery
#     recovery_intensity=0 #recovery is dark
#     rise_time=1 #10 ms for the light to rise
#     time_units='seconds' 
#     point_frequency=100 #start with a frequency of 1000 points per subtrace
#     repeat_cycles=1 #do this once
#     wave=generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity,
#                         pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
#                         rise_time, time_units, point_frequency, repeat_cycles)
#     light_pattern['single_square_20_min_0_max']=wave
    
#     return(light_pattern)

"""
optimized_time_split splits the simulation into small snippets, each with a constant
light intensity, which are simulated in series. This is needed to prevent the sim-
ulations from becoming too stiff for odeint.

"""
def optimized_time_split(test_times_and_light, max_light_change, points_per_segment):   
    test_times_array=test_times_and_light[0]
    test_light=test_times_and_light[1]
    split_points=[] 
    split_begin=0
    sub_arrays_time=[] 
    sub_arrays_light=[]
    split_points.append(0) #start with the zero index split_points

    for i in range(1,len(test_times_array)+1): #test_times_array contains the waveform that will be transformed into
                                                #the appropriate set of sub traces
                                                #the loop progressively increased the length of the selected
                                                #region until the change in light intensity is above the 
                                                #threhold set by max_light_change
        test_range=test_light[split_begin:i]
        if np.max(test_range)- np.min(test_range) > max_light_change: #we have reached the point where the 
                                                                        #light change is larger than the max
            split_points.append(i)
            split_begin=i
    if split_begin<len(test_times_array):
        split_points.append(len(test_times_array)-1)
    for ii in range(1,len(split_points)):
        ptre=split_points[ii]
        ptrb=split_points[ii-1]
        temp_x=np.linspace(test_times_array[ptrb], test_times_array[ptre], points_per_segment)
        #average_light=np.mean(test_light[ptrb:ptre])    #at first, I used the average light intensity over the 
                                                        #entire subtrace, but this was a bad idea because if the 
                                                        #trace was short, it could result in setting the dark baseline
                                                        #to something above zero!
        beginning_light=test_light[ptrb]
        use_this_light=beginning_light
        temp_y=np.linspace(use_this_light, use_this_light, points_per_segment)
        sub_arrays_time.append(temp_x)
        sub_arrays_light.append(temp_y)
    #print('sub_arrays, split_points = ' + str(len(sub_arrays_light)) + ' ' + str(len(split_points)))
    return([sub_arrays_time, sub_arrays_light]) #, split_points])

#def print_constants_table(K):
#    print("{:<30} {:<10}".format('Parameter','Value'))
#    for v in K.items():
#        label, num = v
#        print( "{:<30} {:<10}".format(label, num))

def make_variable_light_constants_set_and_trace_times(K, sub_arrays_time_and_light):
    #K.light_per_L=22
    #print(K.light_per_L)

    sub_arrays_time=sub_arrays_time_and_light[0]
    sub_arrays_light=sub_arrays_time_and_light[1]
    trace_times=[]
    constants_set=[]
    for i in range(len(sub_arrays_time)):
        #print('h ' + str(i))
        K.light_per_L=sub_arrays_light[i][0]
        constants_set.append(K.as_tuple())
        duration=sub_arrays_time[i][-1]-sub_arrays_time[i][0]
        number_of_steps=len(sub_arrays_time[i])
        trace_times.append(np.linspace(0, duration, number_of_steps)) 

    #print('there are ' + str(len(constants_set)) + ' subsets in this trace')
    return([constants_set, trace_times])

def generate_sin_wave(total_duration, max_PAR, light_frequency, points_per_second):
    test_number_points=total_duration*points_per_second
    times_array=np.linspace(0, total_duration, test_number_points, dtype=float)
    sin_wave=[]
    #print('length of test array is: ' + str(len(times_array)))
    for i in times_array:
        sinLt=np.sin(i*2*light_frequency*np.pi-(np.pi/2))
        #sinLt=sinLt+(PAR_offset/max_PAR) #add the offset value, as a fraction of the max_PAR
        sin_wave.append(sinLt)
    sin_wave=np.array(sin_wave)
    sin_wave=sin_wave-np.min(sin_wave)
    if np.max(sin_wave)>0:
        sin_wave=sin_wave/np.max(sin_wave)
    sin_wave=sin_wave*max_PAR
    return([times_array, sin_wave])
