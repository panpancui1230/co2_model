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

def generate_sin_based_light_sequence (total_duration, time_units, max_PAR, 
                                    wave_form, light_frequency, 
                                    point_frequency):
                                    
    #number_segments is the number of segments to split the sequence. The time_slice 
    #will be adjusted for each segment to keep the changes under a certain value.  
    #total_duration, time_units, time_slice, max_PAR, PAR_offset, clipping, light_frequency

    if time_units=='seconds':
        time_div=1
    elif time_units=='minutes':
        time_div=60
    elif time_units=='hours':
        time_div=60*60

    if wave_form=='sin':
        clipping=[0,2]
        PAR_offset=1000

    elif wave_form=='square':
        clipping=[-0.1,0.1] #the numbers used here define the sharpness of the 'square wave'
        PAR_offset=0

    total_duration_in_seconds=total_duration*time_div
    test_number_points=total_duration_in_seconds*point_frequency
    test_times_array=np.linspace(0, total_duration_in_seconds, test_number_points, dtype=float)

    #make the full waveform at high resolution
    test_sin_light_list=[]
    #print('length of test array is: ' + str(len(test_times_array)))
    for i in test_times_array:
        sinLt=np.sin(i*2*np.pi*light_frequency-(np.pi/2))
        sinLt=sinLt+(PAR_offset/max_PAR) #add the offset value, as a fraction of the max_PAR
        if sinLt<clipping[0]: #cannot have negative light, so when the curve goes below, assume it is night
            sinLt=(PAR_offset/max_PAR)-1
        if sinLt>clipping[1]: #gives sharp, square-wave-like limit #cannot have light, so when the
                            #curve goes below, assume it is night
            sinLt=1+(PAR_offset/max_PAR)
        test_sin_light_list.append(sinLt)
        
    test_sin_light_list=np.array(test_sin_light_list)
    test_sin_light_list=test_sin_light_list/np.max(test_sin_light_list)
    test_sin_light_list=test_sin_light_list*max_PAR
    #print(len(test_sin_light_list))
    
    return([test_times_array, test_sin_light_list])


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




#generate a light sequence that contains fluctuations
                    
def fluctuating(total_time, frequency_of_fluctuations, rise_time, max_light, envelope):
    
    #random fluctuations
    #general terms for the pulsed wave:
    duration_of_fluctuation=1/frequency_of_fluctuations
    number_of_cycles=int(total_time/ duration_of_fluctuation)
    
    #set up the basic structure of the square pulse 
    baseline_duration=0.25*duration_of_fluctuation #in seconds
    pulse_duration=0.5*duration_of_fluctuation #100 seconds pulse
    recovery_duration = 0.25*duration_of_fluctuation #100 seconds recovery
    rise_time=1 #1 s for the light to rise
    time_units='seconds' 
    point_frequency=100 #start with a frequency of 1000 points per subtrace
    repeat_cycles=1 #do this once
    rx=np.array([])
    ry=np.array([])
    for i in range(0,number_of_cycles):
        if (envelope == 'sin'):
            sinLt=max_light*(np.sin(i*2*np.pi/number_of_cycles-(np.pi/2))+1)
        else:
            sinLt=max_light
            
        inten=np.random.uniform(0, sinLt, size=3)
        pulse_intensity=inten[0] #pulse is 300 uE m-2 s-1 units
        baseline_intensity=inten[1] #dark baseline
        recovery_intensity=inten[2] #recovery is dark
    
    
        xy=generate_square_wave_based_light_sequence(baseline_duration, baseline_intensity,
                            pulse_duration, pulse_intensity, recovery_duration, recovery_intensity, 
                            rise_time, time_units, point_frequency, repeat_cycles)
                            
        if len(rx)>0:
            rx=np.append(rx, rx[-1]+np.array(xy[0]))
        else:
            rx=np.append(rx, np.array(xy[0]))
        ry=np.append(ry, np.array(xy[1]))
    return [rx,ry]


                    
def generate_square_wave_based_light_sequence (baseline_duration, baseline_intensity, pulse_duration, pulse_intensity, 
                                        recovery_duration, recovery_intensity, rise_time, time_units, point_frequency,
                                        repeat_cycles):
    pulse_times=[] #array that contains the complete time sequence
    pulse_light=[] #array to contain the complete intensity sequence 
    if time_units=='seconds':
        time_div=1
    elif time_units=='minutes':
        time_div=60
    elif time_units=='hours':
        time_div=60*60

    baseline_duration=baseline_duration*time_div
    baseline_points=baseline_duration*point_frequency #calculate the number of points in the baseline
    baseline_time=np.linspace(0, baseline_duration, baseline_points) #generate the baseline array, starting at zero
    baseline_intensity_array= np.linspace(baseline_intensity, baseline_intensity,baseline_points) #fill baseline array with baseline light intensity

    pulse_times=np.append(pulse_times, baseline_time) 
    pulse_light=np.append(pulse_light, baseline_intensity_array)

    riser_duration=rise_time*time_div
    riser_points=riser_duration*point_frequency
    riser_start_time = (baseline_points+1) / point_frequency
    
    riser_end_time = riser_start_time + riser_duration
    riser_time=np.linspace(riser_start_time, riser_end_time, riser_points)
    riser_light=np.linspace(baseline_intensity, pulse_intensity, riser_points)
    #print('rst= ' + str(riser_start_time))
    #print('ret= ' + str(riser_end_time))
    
    pulse_times=np.append(pulse_times, riser_time)
    pulse_light=np.append(pulse_light, riser_light)
    
    pulse_duration=pulse_duration*time_div
    pulse_points=pulse_duration*point_frequency
    #pulse_start_time = (baseline_points + riser_points +1)/point_frequency
    
    pulse_start_time = pulse_times[-1] + 1/point_frequency
    
    pulse_end_time = pulse_start_time + pulse_duration
    
    #print('pst= ' + str(pulse_start_time))
    #print('pet= ' + str(pulse_end_time))


    pulse_time=np.linspace(pulse_start_time, pulse_end_time, pulse_points)
    pulse_light_array=np.linspace(pulse_intensity, pulse_intensity, pulse_points)
    pulse_times=np.append(pulse_times, pulse_time)
    pulse_light=np.append(pulse_light, pulse_light_array)
    
    falling_duration=rise_time*time_div
    falling_points=riser_duration*point_frequency
    
    #falling_start_time = (baseline_points + riser_points + pulse_points + 1) / point_frequency
    falling_start_time = pulse_times[-1] + 1/point_frequency

    falling_end_time = falling_start_time + falling_duration
    
    #print('fst= ' + str(falling_start_time))
    #print('fet= ' + str(falling_end_time))


    falling_time=np.linspace(falling_start_time, falling_end_time, falling_points)
    falling_light=np.linspace(pulse_intensity, recovery_intensity, falling_points)
    
    pulse_times=np.append(pulse_times, falling_time)
    pulse_light=np.append(pulse_light, falling_light)
    
    recovery_duration=recovery_duration*time_div
    recovery_points=recovery_duration*point_frequency
    recovery_start_time = pulse_times[-1] + 1/point_frequency
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
    
def flat_light(day_length_h, PAR, points_per_second):
    day_length_s=day_length_h*60*60
    time_axis_s=np.linspace(0,day_length_s, day_length_s*points_per_second)
    flat_envelope=np.linspace(PAR,PAR, day_length_s*points_per_second)
    return(time_axis_s, flat_envelope)

def sin_light(day_length, max_PAR, points_per_second): 
    
    day_length_s=day_length*60*60
    time_axis_s = np.linspace(0, day_length_s, day_length_s*points_per_second, endpoint=False)
    time_axis_h=time_axis_s/(60*60)
    
    #generate the envelope trace
    envelope=signal.cosine(day_length_s*points_per_second)*max_PAR
    return([time_axis_s, envelope])


def generate_light_profile(input_envelope, fluctuations): 
    time_axis_s=input_envelope[0]
    day_length_s=time_axis_s[-1]
    points_per_second=len(time_axis_s)/day_length_s
    time_axis_h=time_axis_s/(60*60)
    envelope=input_envelope[1]
    
    light_fluct=np.array([])
    light_array=np.array([])
    if fluctuations['type']=='square':

        time_index=0
        if fluctuations['distribution']=='random':
            #print('random')
            if fluctuations['begin'] == 'beginning':
                start_point=0
            else:
                start_point=int(float(fluctuations['begin'])*3600)*points_per_second  #convert to seconds

            if fluctuations['end'] == 'end':
                end_point=int(len(time_axis_s))
            else:
                end_point=int(float(fluctuations['end'])*3600)*points_per_second #convert to seconds

            while time_index<len(time_axis_s): #day_length_s:
                duration_of_fluctuation=np.random.randint(fluctuations['tao'])
                if (time_index>start_point-1) and (time_index<end_point+1):
                    fluctuation_amplitude=np.random.uniform(float(fluctuations['amplitude'][0]), float(fluctuations['amplitude'][1]))
                    light_value=envelope[len(light_fluct)]*(1-fluctuation_amplitude)
                    light_fluct=np.append(light_fluct, np.linspace(fluctuation_amplitude, fluctuation_amplitude, (duration_of_fluctuation*points_per_second)))
                    light_array=np.append(light_array, np.linspace(light_value, light_value, (duration_of_fluctuation*points_per_second)))

                else:
                    light_fluct=np.append(light_fluct, np.linspace(0, 0, (duration_of_fluctuation*points_per_second)))
                    light_array=np.append(light_array, envelope[time_index:time_index+int(duration_of_fluctuation*points_per_second)])

                    #light_value = envelope[len(light_fluct)]


                time_index=time_index+int(duration_of_fluctuation*points_per_second)
            light_array_smoothed=smooth(light_array, int(fluctuations['smooth_points']))
            out_put_light_array=light_array_smoothed[0:len(time_axis_s)]
    else:
        out_put_light_array=envelope
        
    return([time_axis_s, out_put_light_array])

def multiply_light_pattern(pattern, times_to_repeat):
    wave_x=np.array([])
    wave_y=np.array([])
    time_offset=0
    for index in range(0,times_to_repeat):
        wave_x=np.append(wave_x, pattern[0]+time_offset)
        wave_y=np.append(wave_y, pattern[1])
        time_offset=wave_x[-1]
    return [wave_x, wave_y]
    
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

    

#*******************************************************************************
#*******************************************************************************
#                   Code related to ODE simulations                            *
#*******************************************************************************
#*******************************************************************************


"""

Notes on the functions calc_K_b6f and calc_v_b6f:
    To estimate the effective rate of PQH2 oxidation at the cytochrome b6f complex, 
    we need to consider the redox states of PQH2, PC as well as the Dy and DpH. 
    Because of the Q-cycle, 2 H+ are transferred into the lumen for each electron 
    passed from PQH2 to PC, i.e. 
    
        0.5PQH2 + b6f(protonated) + PC(ox) --k_b6f--> PQ + b6f(protonated) + PC(red) + 2Hin
    
    The forward rate constant is k_b6f, but the reaction is reversible, so that
    
        0.5PQH2 + b6f(protonated) + PC(ox) <--k_b6f_reverse-- PQ + b6f(protonated) + PC(red) + 2Hin
    
    k_b6f_reverse is a function of pmf because the Q-cycle in the forward direction works against 
    both DpH and Dy. (Note that this thermodynamic effect is in addition to the kinetic effect on 
    the deprotonation of the Rieske protein.) We simplify this for the simulation as follows:
    
        Keq_b6f = Em(Pc_ox/PC_red) - Em(PQ/PQH2) - pmf

    In other words, the eqiulibirum constant is determined by the redox potentials of the donor 
    and acceptor together and the pmf. We use unity as the scaling factor for the pmf 
    contributions becaus ecause one proton translocated to the lumen per e- t
    ransferred (together with one e- charge moved from the p- to the n-side) equilibrium.
    
        k_b6f_reverse = k_b6f / Keq
    
    In principle we could simulate the effects of changing PQH2 and PC redox states in two ways, 
    either using the simulated concentrations of PQH2 and PC together with the standard E'0 values, 
    or accounting for the concentrations in the Em values. We chose the former because 
    it better fits the form of the ODE equations and is a bit simpler to calculate. Thus,
    
        v_b6f=[PQH2][PC_ox]k_b6f - [PQ][PC_red]k_b6f_reverse

    where E'0(Pc_ox/PC_red) = 0.370 V, pH-independent under our conditions; E'0(PQ/PQH2) = 0.11 V at pH=7, 
    but pH-dependent so that: 
        
        E'0(PQ/PQH2) = 0.11 V - (7-pHlumen) * 0.06
        
    at pH=7, but pH-dependent so that:

        Keq_b6f = E'0(Pc_ox/PC_red) - E'0(PQ/PQH2) - pmf = 0.370 - 0.11 + .06 * (pHlumen - 7.0) - pmf

    So, the full set of equations is:
        Em7_PC=0.37 Em_PC=Em7_PC Em7_PQH2 = 0.11 Em_PQH2= Em7_PQH2 + 0.06*(pHlumen-7.0)
    Keq_b6f = 10**((Em_PC - Em_PQH2 - pmf)/.06)
    k_b6f_reverse = k_b6f / Keq
    v_b6f=PQH2PC_oxk_b6f - PQPC_redk_b6f_reverse

"""


def calc_k_b6f(max_b6f, b6f_content, pHlumen, pKreg):    
    #pHmod is the fraction of b6f complex that is deprotonated
    pHmod=(1 - (1 / (10 ** (pHlumen - pKreg) + 1)))
    b6f_deprot=pHmod*b6f_content
    k_b6f=b6f_deprot * max_b6f
    return(k_b6f)

#v_b6f=calc_v_b6f(max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf)

def calc_v_b6f(max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf):    
    pHmod=(1 - (1 / (10 ** (pHlumen - pKreg) + 1)))
    b6f_deprot=pHmod*b6f_content

    Em_PC=Em7_PC
    Em_PQH2= Em7_PQH2 - 0.06*(pHlumen-7.0)

    Keq_b6f = 10**((Em_PC - Em_PQH2 - pmf)/.06)
    k_b6f=b6f_deprot * max_b6f 

    k_b6f_reverse = k_b6f / Keq_b6f
    #print('Keq for PQH2 to PC + pmf is: ' + str(Keq_b6f))
    f_PQH2=PQH2/(PQH2+PQ) #want to keep the rates in terms of fraction of PQHs, not total number
    f_PQ=1-f_PQH2
    v_b6f=f_PQH2*PC_ox*k_b6f - f_PQ*PC_red*k_b6f_reverse 
    return(v_b6f)
"""
Notes on the NDH activity for CEF
  Forward reaction
	Fd_red + 0.5 PQ + 3H_stroma --k_NDH--> Fd_ox + 0.5 PQH2 + 2H_lumen

  Reverse reaction
	Fd_red + 0.5 PQ + 3H_stroma <--k_NDH_reverse-- Fd_ox + 0.5 PQH2 + 2H_lumen

Em_Fd = -0.42
Em_PQH2_7 = 0.11
Em_PQH2 = 0.11 - 0.06*(pHstroma-7.0), # T = 20 degree C

deltaEm = Em_PQH2 - Em_Fd = 0.53 - 0.06*(pHstroma-7.0)

deltaG_NDH = z*F*deltaEm + 2*F*pmf, 
here z = -1, and 2 is for the number of protons pumped by NDH, so
Keq_NDH = 10 **(((0.53-0.06*(pHstroma-7.0)-2*pmf)/0.06)

Assuming k_NADH = 100/s
k_NDH_reverse = k_NDH /Keq_NDH

v_NDH = k_NDH*Fd_red*f_PQ - k_NDH_reverse*Fd_ox*f_PQH2

For each e- transferred, d_charge = 2, dH_lumen = 2, dH_stroma = -3
"""
def calc_v_NDH(Em_Fd, Em7_PQH2, pHstroma, pmf, k_NDH, Fd_red, Fd_ox, PQ, PQH2):
    Em_PQH2 = Em7_PQH2 - 0.06*(pHstroma - 7.0)
    deltaEm = Em_PQH2 - Em_Fd
    Keq_NDH = 10**((deltaEm - pmf*2)/0.06)
    k_NDH_reverse = k_NDH/Keq_NDH
    #f_PQ = PQ/(PQ+PQH2)
    #f_PQH2 = 1.0-f_PQ
    v_NDH = k_NDH*Fd_red*PQ - k_NDH_reverse*Fd_ox*PQH2
    return (v_NDH)
def calc_v_PGR(PGR_vmax, Fd_red, PQ, PQH2):
    #Fd_red + 1/2PQ + H+_Stroma--> Fd_ox +1/2 PQH2, Hill coefficient for PGR assumed at 4
    #without considering back reaction,this is the guess of PGR5/PGRL1 cyclic etr
    v_PGR = PGR_vmax * (Fd_red**4/(Fd_red**4+0.1**4))*PQ/(PQ+PQH2)
    #The reason 0.1 was chose is that Fd_red level does not seem to go over 0.2
    return v_PGR
#calculate the rate of V<-- -->Z reactions, assuming a pH-dependent VDE and a pH-independent ZE
def calc_v_VDE(VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pHlumen, V, Z):    

    #VDE_Hill is the Hill coefficient for the VDE reaction
    #pKvde is the pKa for protonation of VDE
    #VDE_max_turnover_number is the maximum turnover rate (at low pH for VDE)
    #kZE is the rate constant for the ZE reaction

    #pHmod is the fraction of VDE complex that is deprotonated
    pHmod= 1 / (10 ** (VDE_Hill*(pHlumen - pKvde)) + 1)
    
    #pHmod=1-(1 - (1 / (10 ** (VDE_Hill*(pHlumen - pKvde) + 1))))
    #pHmod=(1-(1 / (10 ** (VDE_Hill*(pHlumen - pKvde) + 1))))
    #print(pHmod)
    #calculate the net change in Z
    v_Z = V* VDE_max_turnover_number*pHmod - Z*kZE
    v_V = -1* v_Z
    
    return(v_Z, v_V)
    
#calculate the rate of V<-- -->Z reactions, assuming a pH-dependent VDE and a pH-independent ZE
def calc_PsbS_Protonation(pKPsbS, pHlumen):    

    #VDE_Hill is the Hill coefficient for the VDE reaction
    #pKvde is the pKa for protonation of VDE
    #VDE_max_turnover_number is the maximum turnover rate (at low pH for VDE)
    #kZE is the rate constant for the ZE reaction

    #pHmod is the fraction of VDE complex that is deprotonated
    PsbS_H=1 / (10 ** (3*(pHlumen - pKPsbS)) + 1)
    
    return(PsbS_H)

"""
However, one may also consider that there is a maximal (saturating turover rate 
(saturation point), as shown by Junesch and Grabber (1991)
http://dx.doi.org/10.1016/0014-5793(91)81447-G
Their data shows a roughly n=1 pmf-dependence, similar to a pH titration curve, but for pmf, which can 
be simulated by changing the code to include this term.
    
"""
def ATP_synthase_actvt(t):#based on gH+ data
    x = t/T_ATP
    actvt = 0.2 + 0.8*(x**4/(x**4 + 1))
    #x = t/174.5
    #nth = (x-1)*3.1
    #actvt = 0.112 + 0.888*1/(1+np.exp(-nth))
    #actvt can be treated as pmf responsive fraction and 1-actvt is pmf_inert fraction
    #this actvt may due to ATP synthase oxidized-->reduced delay or/and ATP/ADP, Pi limitation

    return actvt
    
def Vproton_pmf_actvt(pmf, actvt, ATP_synthase_max_turnover, n):# fraction of activity based on pmf, pmf_act is the half_max actvt pmf
    v_proton_active = 1 - (1 / (10 ** ((pmf - 0.132)*1.5/0.06) + 1))#reduced ATP synthase
    v_proton_inert = 1-(1 / (10 ** ((pmf - 0.204)*1.5/0.06) + 1))#oxidized ATP synthase
    
    v_active = actvt * v_proton_active * n * ATP_synthase_max_turnover
    v_inert = (1-actvt) * v_proton_inert * n * ATP_synthase_max_turnover
    
    v_proton_ATP = v_active + v_inert

    #the factor 1.5 is used as a hill coefficient for ATP synthase, adjusted from reference Fig.3
    #Note that the above experiments were done at deltaGatp = 30kJ/mol + RTln(1.2/5)
    #so pmf_act needs to be adjusted by the following function ATP_deltaG into pmf_addition
    return (v_proton_ATP)
    #the following is another way to simulate, but it cannot realize observed data
    #Patp_red = 1/(1+np.exp(-t+6))
    #Patp_ox = 1- Patp_red
    #pmf_act_red = 0.132 #see the following reference, dpH_half = 2.2 and 3.4 for red and ox
    #pmf_act_ox = 0.204#Ulrike Junesch and Peter Graber BBA 893(1987) 275-288
    #pmf_addition = ATP_pmf_addition(ATP_pool, ADP_pool, Pi)
    #pmf_addition is the deltaG difference between realtime deltaGatp and reference
    #experimental deltaGatp, and then converted to pmf units
    #pmf_act_red = pmf_act_red + pmf_addition
def V_H_light(light_per_L, v_proton_ATPase, pmf, Hlumen, k_leak = 3*10**7):
    if light_per_L>0.0:
        V_H = v_proton_ATPase + pmf*k_leak*Hlumen
    else:
        V_H = pmf*k_leak*Hlumen# this term is used for dark relaxation,
        #ATP synthase actvt dependent but does not make ATP
    return V_H

def calc_CBC_NADPH(k_CBC,t,d_ATP_made):
    NADPH_CBC_t =  k_CBC*(1.0-np.exp(-t/600))
    NADPH_CBC_ATP =0.6*d_ATP_made
    NADPH_CBC = min([NADPH_CBC_ATP,NADPH_CBC_t])
    return NADPH_CBC
        

def Calc_Dy(CIONlumen,CIONstroma, AIONlumen,AIONstroma):
    Dy_ION = 0.06*np.log10((CIONlumen+0.5*AIONstroma)/(CIONstroma+0.5*AIONlumen))
    #absolute value of electrical potential, the 0.5 is due to the PCl- = 0.5 PK+
    return Dy_ION

def calc_pmf_act(ATP_pool, ADP_pool, Pi):
    DeltaGatp_zero =  36.0#Petersen et al. 2012 PNAS, comparsion of the H+/ATP ratios of mF0F1 cF0F1 ATP synthase
    DeltaDeltaGatp =2.44 * np.log(ATP_pool/(ADP_pool*Pi))#np.log is natural log
    #the pH effect in stroma is ignored at this stage
    DeltaGatp_KJ_per_mol = DeltaGatp_zero + DeltaDeltaGatp
    #"Taras K. Antal • Ilya B. Kovalenko, Andrew B. Rubin • Esa Tyystjarvi, Photosynthesis-related quantities for education and modeling. Photosynth Res (2013) 117:1–30
    #D. Heineke et al., Redox transfer across the inner chloroplast envelope membrane. Plant Physiol. 95, 1131–1137 
    #(1991). DeltaGatp_KJ_per_mol=50.0, this number should be under the light adapted conditions"
    #Assuming Pi is 1.5 mM constant, under light ATP is 1 mM, then the ADP under light is 0.184 mM
    #This agrees with report of ATP/ADP about 5 under the light adapted condition
    #This gives a pool of ATP of 7/PSII, ADP 1.3/PSII under the light
    #ATP/ADP = 1 for dark adapted plant.
    
    #convert DGATP into volts
    pmf_act = DeltaGatp_KJ_per_mol/(96.485*4.667)#Faraday constant
    #pmf_addition = DeltaDeltaGatp/(96.485*4.667) #Faraday constant
    return (pmf_act)

"""
# Calc_Phi2 gives an estiamte of Phi2 based on QA redox state
# and NPQ. The dertivation of the equation is based on that described in


D.M. Kramer, G. Johnson, O. Kiirats, G.E. Edwards (2004) New fluorescence 
parameters for the determination of QA redox state and excitation energy fluxes. 
Photosynthesis research 79, 209-218.

and

S. Tietz, C.C. Hall, J.A. Cruz, D.M. Kramer (2017) NPQ(T): a chlorophyll fluorescence p
arameter for rapid estimation and imaging of non-photochemical quenching of excitons i
n photosystem II associated antenna complexes. Plant Cell and Environment In Press.


The following is a derivation for determining Phi2 from NPQ and QA redox state
Recall that NPQ is the ratio of:
    
    NPQ=kNPQ/(kd+kf)

    and

    kNPQ=NPQ/(kf+kd)

    where kNPQ if the rate constant for NPQ, kq is the k for intrinsic non-radiative, and kf is the rate constant for fluorescence

The maximal PSII quantum yield is:

    Phi2_max = kpc/(kd + kf +kpc) = 0.83

    where kpc is the maximal rate contant for quenching of excitation energy by open PSII
and thus: 
    Phi2 = QAkpc/(kf + kd + kNPQ + QAkpc) = QAkpc/(kf + kd + NPQ/(kf+kd) + QAkpc) 
            = QAkpc/(kf + kd + NPQ/(kf+kd) + QAkpc)
        
    1/Phi2= (kd + kf + kNPQ + QAkpc)/QAkpc = 1 + (kd+kf+ kNPQ)/QAkpc
            = 1 + (kf+kd)/QAkpc + kNPQ/QAkpc 1/(PHI2(kf+kd)) 
            = 1/(kf+kd) + 1/(QAkpc) + kNPQ/((kf+kd)QAkpc) 1/(PHI2(kf+kd)) 
            = 1/(kf+kd) + 1/(QAkpc) + NPQ/QAkpc
            = 1/(kf+kd) + 1/(QA*kpc) + NPQ/(QA*kpc)
            
    1/Phi2_max = (kd + kf)/kd +kpc/kpc 
                = 1+ (kf + kd)/kpc 0.83-1 
                = (kf + kd)/kpc =0.17 kpc/(kf+kd)=5.88
                
    kpc/(PHI2(kf+kd)) = kpc/(kf+kd) + kpc/(QAkpc) + kpcNPQ/QAkpc
    
    5.88/Phi2 = 5.88 + 1/QA + NPQ/QA = 5.88 + (1+NPQ)/QA 1/Phi2=1 + (1+NPQ)/(5.88*QA)
    Phi2=1/(1+(1+NPQ)/(5.88*QA))

   = 1
   ____________
   1+ (1+NPQ)
      _______
      5.88*QA


"""

def Calc_Phi2(QA, NPQ):
    Phi2=1/(1+(1+NPQ)/(4.88*QA))
    return Phi2

"""
Calc_PhiNO_PhiNPQ gives estiamtes of PhiNO and PhiNPQ based on the
equations is based on that described in

D.M. Kramer, G. Johnson, O. Kiirats, G.E. Edwards (2004) New fluorescence 
parameters for the determination of QA redox state and excitation energy fluxes. 
Photosynthesis research 79, 209-218.

S. Tietz, C.C. Hall, J.A. Cruz, D.M. Kramer (2017) NPQ(T): a chlorophyll fluorescence 
parameter for rapid estimation and imaging of non-photochemical quenching of excitons in 
photosystem II associated antenna complexes. Plant Cell and Environment In Press.

and derived using the approach detailed for Phi2.

"""

def Calc_PhiNO_PhiNPQ(Phi2, QA, NPQ):
    PhiNO=1/(1+NPQ + ((Phi2+NPQ)/(1-Phi2)))
    PhiNPQ=1-(Phi2+PhiNO)
    return PhiNO, PhiNPQ


"""
Notes on calculation of PSII recombination rates:
    
    I used the equations presented in Davis et al. 2016 

G.A. Davis, A. Kanazawa, M.A. Schöttler, K. Kohzuma, J.E. Froehlich, A.W. 
Rutherford,M. Satoh-Cruz, D. Minhas, S. Tietz, A. Dhingra, D.M. Kramer 
(2016) Limitations to photosynthesis by proton motive force-Induced 
photosystem II photodamage eLife eLife 2016;5:e16921.

Specifically,there are two parts to the estimation of rates of recombination:
    
    v_recombination = k_recomb*QAm*(10**(Dy/.06) + fraction_pH_effect*10**(7.0-pHlumen))

    where k_recomb is the rate constant for S2QA- recombination in the 
    absence of a field in s^-1. Dy is the delta_psi in volts, QAm is the 
    content of reduced QA~0.33 or one recombination per 3 s, s seen in the 
    presence of DCMU. The term fraction_pH_effect rtepresents the fraction 
    of S-states that are both pH-sensitive i.e. involve release of protons, 
    and are unstable (can recombine).

Then, 10**(7.0-pHlumen) represents the change in equilibrium constant
for the Sn+1 P680 <--> SnP680+, as a result of changes in lumen pH 

"""

def recombination_with_pH_effects(k_recomb, QAm, Dy, pHlumen, fraction_pH_effect):
    delta_delta_g_recomb= Dy + .06*(7.0-pHlumen)
    v_recomb = k_recomb*QAm*10**(delta_delta_g_recomb/.06)
    
    #v_recomb = k_recomb*QAm*(10**((Dy/.06) + fraction_pH_effect*10**(7.0-pHlumen)))        
    return(v_recomb)

def Cl_flux_relative(v):
    Cl_flux_v = 332*(v**3) + 30.8*(v**2) + 3.6*v
    #relative to Cl flux thru VCCN1. when driving force is 0.1 Volt,
    #Cl_flux_v is 1. empirical equation was obtained from
    # Herdean et al. 2016 DOI: 10.1038/ncomms11654
    return Cl_flux_v
def KEA_reg(pHlumen, QAm):
    qL = 1-QAm
    qL_act = qL**3/(qL**3+0.15**3)
    pH_act =1/(10**(1*(pHlumen-6.0))+1)
    f_KEA_act = qL_act * pH_act
    return f_KEA_act
#Function f calculates the changes in state for the entire systems
def square_light(t, light_intensity, duration = 20, unit= 'min',\
                 t0 = 0, prior_light = 0, post_light = 0 ):
    """
    caculate a light intensity in a square wave

    Parameters
    ----------
    t : int or float
        time point in seconds.
    light_intensity : int or float, unit in umol photons/m^2/s
        light intenity during given duration: t0 to t0+duration.
    duration : numeric, optional
        DESCRIPTION. The default is 20.
    unit : string, optional
        DESCRIPTION. The default is 'min'.
    t0 : numeric, optional
        DESCRIPTION. The default is 0 (seconds).
    prior_light : numeric, optional
        LIGHT INTENSITY BEFORE t0. The default is 0.
    post_light : numberic, optional
        LIGHT INTENSITY AFTER DURATION. The default is 0.

    Returns
    -------
    par : light intensity of a given time point

    """
    if unit in ['m', 'minutes', 'min', 'minute']:
        duration = duration *60#convert to seconds
    if unit in ['hours','h','hour', 'hr', 'hrs']:
        duration = duration *3600#convert to seconds
    if t < t0:
        par = prior_light
    elif t < (t0+duration):
        par = light_intensity
    else:
        par = post_light
    return par


def sin_light_fluctuating(t, freq, PAR_max, PAR_min, t0=0, PAR_0=None ):
    """
    calculate light intensity when given t

    Parameters
    ----------
    t : numeric, float, or np.array of numeric
        TIME POINTS in seconds.
    freq : float
        FREQUENCY of sine wave.
    PAR_max : numeric
        Maximum PAR.
    PAR_min : float
        Minimum PAR.
    t0 : float, optional
        initial time. The default is 0.
    PAR_0 : None or float, optional
        initial PAR. The default is None.

    Returns
    -------
    par: light intensity of at a time point.

    """

    if PAR_0 == None:
        PAR_0 = (PAR_max+PAR_min)/2
    elif PAR_0 < PAR_min or PAR_0 > PAR_max:
        PAR_0 = (PAR_max+PAR_min)/2
        #print('The given PAR_0 is out of range, (PAR_max+PAR_min)/2 was used')
    A = (PAR_max - PAR_min)/2#amplitude of light fluctuation
    sin_phi = (PAR_0-(PAR_max+PAR_min)/2)/A# sin(2π*freq*0+phi)
    phi = np.arcsin(sin_phi)#initial phase
    par = PAR_min + A * (1+ np.sin(2*np.pi*freq*(t-t0) + phi))
    return par      

def light(t, Duration_T0, par0, frequency, par_max, par_min):
    if t <= Duration_T0:
        par = square_light(t, par0, Duration_T0, 'seconds')
    else:
        par = sin_light_fluctuating(t, frequency, par_max, par_min, \
                                    t0 = Duration_T0, PAR_0 = par0)
    return par


def f(t, y, pKreg, max_PSII, kQA, max_b6f, lumen_protons_per_turnover, PAR, ATP_synthase_max_turnover, 
    PSII_antenna_size, Volts_per_charge, perm_K, n, Em7_PQH2, Em7_PC,Em_Fd, PSI_antenna_size, 
    buffering_capacity, VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pKPsbS, max_NPQ, k_recomb, k_PC_to_P700, 
    triplet_yield, triplet_to_singletO2_yield, fraction_pH_effect, k_Fd_to_NADP, k_CBC, k_KEA, k_VCCN1, k_CLCE, k_NDH): 
    
    #The following are holders for paramters for testing internal functions of f
    PAR = light(t, 1200, LIGHT, FREQUENCY, 900, 100)
    light_per_L=0.84 * PAR/0.7
    #we got 4.1 nmol Chl per 19.6 mm^2 leaf disc, which translate into 210 umol Chl/m2
    #210 umol Chl/m2, PSII/300 Chl ==> 0.7 umol PSII/m2, ==>(PAR/0.7) photons/PSII
    ###So the light_per_L means the photons/PSII that hit all the thylakoid membranes, and absorbed by the leaf
    #the 0.84 is the fraction of light a leaf typically absorbs

    
    QA, QAm, PQ, PQH2, Hin, pHlumen, Dy, pmf, deltaGatp, Klumen, Kstroma, ATP_made,\
    PC_ox, PC_red, P700_ox, P700_red, Z, V, NPQ, singletO2, Phi2, LEF, Fd_ox, Fd_red,\
    ATP_pool, ADP_pool, NADPH_pool, NADP_pool,Cl_lumen, Cl_stroma, Hstroma, pHstroma =y
    
    PSII_recombination_v=recombination_with_pH_effects(k_recomb, QAm, Dy, pHlumen, fraction_pH_effect)
        
    dsingletO2=PSII_recombination_v*triplet_yield*triplet_to_singletO2_yield

    #calculate pmf from Dy and deltapH 
    pmf=Dy + 0.06*(pHstroma-pHlumen)

    #***************************************************************************************
    #PSII reations
    #****************************************************************************************
    #first, calculate Phi2
    Phi2=Calc_Phi2(QA, NPQ) #I use the current' value of NPQ. I then calculate the difference below 

    #calculate the number of charge separations in PSII per second
    PSII_charge_separations=PSII_antenna_size*light_per_L * Phi2
    
    #The equilibrium constant for sharing electrons between QA and the PQ pool
    #This parameter will be placed in the constants set in next revision
    
    Keq_QA_PQ=200
    
    #calculate the changes in QA redox state based on the number of charge separations and equilibration with 
    #the PQ pool
    dQAm = PSII_charge_separations  + PQH2*QA*kQA/Keq_QA_PQ  - QAm * PQ * kQA - PSII_recombination_v
    dQA = -1*dQAm

    #***************************************************************************************
    #PQ pool and the cyt b6f complex
    #***************************************************************************************

    #vb6f = k_b6f(b6f_max_turnover_number, b6f_content, pHlumen, pKreg, PQH2)

    #b6f_content describes the relative (to standaard PSII) content of b6f 
    #This parameter will be placed in the constants set in next revision
    b6f_content=0.433 #Journal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014
    #doi:10.1093/jxb/eru090 Advance Access publication 12 March, 2014
    #Mathias Pribil1, Mathias Labs1 and Dario Leister1,2,* Structure and dynamics of thylakoids in land plantsJournal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014
   
    #calc_v_b6f return the rate of electron flow through the b6f complex
    v_b6f=calc_v_b6f(max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf)
    
    v_NDH = calc_v_NDH(Em_Fd, Em7_PQH2, pHstroma, pmf, k_NDH, Fd_red, Fd_ox, PQ, PQH2)
    d_Hlumen_NDH = v_NDH*2 #change in lumen protons
    d_charge_NDH = d_Hlumen_NDH # change in charges across the membrane
    #d_Hstroma_NDH = v_NDH*3 # change in stroma protons
    
    ##PGR regulation, attempted
    PGR_vmax = 0#It seems this function does not impact the kinetics much.
    v_PGR = calc_v_PGR(PGR_vmax, Fd_red, PQ, PQH2)

    #calculate the change in PQH2 redox state considering the following:
    #PQ + QAm --> PQH2 + QA ; PQH2 + b6f --> PQ    
    PSI_charge_separations= P700_red * light_per_L * PSI_antenna_size * Fd_ox
    #aleternatively,
    #PSI_charge_separations = P700_red*light_per_L*PSI_antenna_size*FB/(FB+FB_minus)
    #PSI_to_Fd = FB_minus*Fd_ox*k_FB_Fd
    #d_FB_minus = PSI_charge_separations-PSI_to_Fd
    #d_FB = -d_FB_minus
    

    dPQH2 = (QAm * PQ * kQA + v_NDH + v_PGR - v_b6f - PQH2*QA*kQA/Keq_QA_PQ)*0.5 
    dPQ = -1*dPQH2

    #***************************************************************************************
    #PSI and PC reactions:
    #***************************************************************************************

    #Calculate the changes in PSI redox state. The current form is greatly simplified, 
    #but does consider the need for oxidized Fd.
    #At this point, I assumed that the total Fd pool is unity
    
    
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
    #alternatively,
    #dFd_red = PSI_to_Fd - k_Fd_to_NADP*Fd_red*NADP_pool - v_NDH-V_me
    
    #***************************************************************************************
    # ATP synthase reactions:
    #***************************************************************************************
    #However, one may also consider that there is a maximal (saturating turover rate 
    #(saturation point), as shown by Junesch and Grabber (1991)
    #http://dx.doi.org/10.1016/0014-5793(91)81447-G
    #    def Vproton(ATP_synthase_max_turnover, n, pmf, pmf_act):
    #    return (ATP_synthase_max_turnover*n*(1 - (1 / (10 ** ((pmf - pmf_act)/.06) + 1))))
    #vHplus=Vproton(ATP_synthase_max_turnover, n, pmf, pmf_act)
    
    #ATP_synthase_driving_force=pmf-(deltaGatp/n) #this is positive if pmf is sufficient to drive 
    #reaction forward, assuming ATP synthase activity is time dependent, derived from gH+ data
    # data courtesy from Geoff and Dave
    #Pi = 0.0025 - ATP_pool/7000
    #pmf_act = calc_pmf_act(ATP_pool, ADP_pool, Pi)
    Hlumen = 10**(-1*pHlumen)
    Hstroma = 10**(-1*pHstroma)
    activity = ATP_synthase_actvt(t)
    
    d_protons_to_ATP = Vproton_pmf_actvt(pmf, activity, ATP_synthase_max_turnover, n)
    d_H_ATP_or_passive = V_H_light(light_per_L, d_protons_to_ATP, pmf, Hlumen)                              
    #d_protons_to_ATP_red = Vproton(ATP_synthase_max_turnover, n, pmf, pmf_act_red)*Patp_red
    #d_protons_to_ATP_ox = Vproton(ATP_synthase_max_turnover, n, pmf, pmf_act_ox)*Patp_ox
    #d_protons_to_ATP = d_protons_to_ATP_red + d_protons_to_ATP_ox
        
    d_ATP_made=d_protons_to_ATP/n                                        
    #The CBC is either limited by Phi2 or by the activation kinetics, take the minimum
    #NADPH_phi_2 = (PSII_charge_separations - PSII_recombination_v)*0.5

    NADPH_CBC = k_CBC*(1.0-np.exp(-t/900))*(np.log(NADPH_pool/NADP_pool)-np.log(1.25))/(np.log(3.5/1.25))#calc_CBC_NADPH(k_CBC, t, d_ATP_made)
    #this number in "np.exp(-t/600)" is important, which impacts the shape of the curves
    dNADPH_pool=0.5 * k_Fd_to_NADP*NADP_pool*Fd_red - NADPH_CBC
    dNADP_pool=-1*dNADPH_pool
    
    dLEF=k_Fd_to_NADP*NADP_pool*Fd_red
    
    d_ATP_consumed = d_ATP_made#NADPH_CBC*5/3 + (ATP_pool/(ADP_pool+ATP_pool)-0.5)*1.2#ATP_pool*(ATP_pool/ADP_pool-1)
    #***************************************************************************************
    #Proton input (from PSII, b6f and PSI) and output (ATP synthase) reactions :
    #***************************************************************************************
    #calculate the contributions to lumen protons from PSII, assuming a average of 
    #one released per S-state transition. In reality, the pattern is not 1:1:1:1, 
    #but for simplicity, I am assuming that the S-states are scrambled under our 
    #illumination conditions. This is described in more detail in the manuscript.
    
    d_protons_from_PSII = PSII_charge_separations - PSII_recombination_v

    #calculate the contributions to Dy from PSII
    charges_from_PSII = PSII_charge_separations - PSII_recombination_v
    
    #calculate the contributions to lumen protons from b6f complex
    #assuming the Q-cycle is engaged, asn thus
    #two protons are released into lumen per electron from
    #PQH2 to PC
    """
    C.A. Sacksteder, A. Kanazawa, M.E. Jacoby, D.M. Kramer (2000) The proton to electron 
    stoichiometry of steady-state photosynthesis in living plants: A proton-pumping Q-cycle 
    is continuously engaged. Proc Natl Acad Sci U S A 97, 14283-14288.

    """
    d_protons_from_b6f = v_b6f*2 #two protons per electron transferred from PQH2 to PC

    #calculate the contributions to Dy from Q-cycle turnover
    #one charge through the low potential b chain per
    #PQH2 oxidized
    charges_from_b6f = v_b6f
     
    #add up the changes in protons delivered to lumen
    #note: net_protons_in is the total number of protons input into the lumen, including both free and bound.
    net_protons_in = d_protons_from_PSII + d_protons_from_b6f + d_Hlumen_NDH - d_H_ATP_or_passive
    #net_protons_stroma = d_protons_to_ATP - v_b6f - d_Hstroma_NDH - QAm * PQ * kQA + PQH2*QA*kQA/Keq_QA_PQ  - dNADPH_pool - d_ATP_made
    #each ATP synthesis consumes one proton

    #see appendix for explanation
    
    #K_deltaG=0.06*np.log10(Kstroma/Klumen) - Dy
    
    #the KEA reaction looks like this:
    # H+(lumen) + K+(stroma) <-- --> H+(stroma) + K+(lumen)
    #and the reaction is electroneutral, 
    #so the forward reaction will depend on DpH and DK+ as:
    
    
    f_actvt = KEA_reg(pHlumen, QAm)
    v_KEA = k_KEA*(Hlumen*Kstroma -  Hstroma*Klumen)*f_actvt#/(10**(2*(pHlumen-6.5))+1)
    
    #Pck = 1/(1+np.exp(39.5*0.66*(-0.003-Dy)))#probability of v_K_channel open
    K_deltaG=-0.06*np.log10(Kstroma/Klumen) + Dy
    v_K_channel = perm_K * K_deltaG*(Klumen+Kstroma)/2
    
   
    #v_K_channel = Pck*perm_K * Dy * 39.5*(Klumen- Kstroma*np.exp(-39.5*Dy))/(1-np.exp(-39.5*Dy))#eq regular
    #v_K_channel = Pck*perm_K * (Klumen*np.exp(39.5*Dy)- Kstroma*np.exp(-39.5*Dy))#eq Hui Lyu
    #Adjusted from Hui Lyu and Dusan Lazar Journal of Theoretical Biology 413 (2017) 11-23, 39.5 = F/RT
    #It seems the flux of K+ is behave similar between Kramer and Lazar simulations.
    #Now the equation considers the  Goldman–Hodgkin–Katz flux equation
    #Hille, Bertil (2001) Ion channels of excitable membranes, 3rd ed.,p. 445, ISBN 978-0-87893-321-1
    
    #Next, use this to calculate a flux, which depends
    #on the permeability of the thylakoid to K+, perm_K:
    net_Klumen =  v_KEA - v_K_channel
    
    #if Dy is +60 mV, then at equilibrium, Kstroma/Klumen should be 10, at which point Keq=1.
    #the Keq is equal to kf/kr, so the rato of fluxes is 

    #net_Klumen=perm_K * K_Keq - perm_K/K_Keq 
    #calculate the change in lumen [K+] by multiplying the change in K+ ions
    #by the factor lumen_protons_per_turnover that relates the standard
    #complex concentration to volume:
    #the term is called "lumen_protons_per_turnover" but is the same for 
    #all species
        
    dKlumen = net_Klumen*lumen_protons_per_turnover
    
    #We assume that the stromal vaolume is large, so there should be 
    #no substantial changes in K+
    
    dKstroma=0
    #########now calculating the movement of Cl- and its impact####

    driving_force_Cl = 0.06* np.log10(Cl_stroma/Cl_lumen) + Dy
    v_VCCN1 = k_VCCN1 * Cl_flux_relative(driving_force_Cl) * (Cl_stroma + Cl_lumen)/2#Cl_flux_relative()
    ##v_VCCN1 is rate of Cl- moving into lumen, v_CLCE is rate of Cl- moving in/out

    #here CLCE is assumed one H+ out, two Cl- comes in
    v_CLCE =  k_CLCE*(driving_force_Cl*2+pmf)*(Cl_stroma + Cl_lumen)*(Hlumen+Hstroma)/4
    #v_CLCE = k_CLCE *(Cl_lumen * Hlumen - Cl_stroma * Hstroma)
    
    net_Cl_lumen_in = v_VCCN1 + 2*v_CLCE
    dCl_lumen = net_Cl_lumen_in * lumen_protons_per_turnover
    dCl_stroma = -0.1*dCl_lumen
    
    #***************************************************************************************
    #Buffering capacity and calculation of lumen pH:
    #***************************************************************************************
    #H_leak = Per_H * ([Hlumen]-[Hstroma])
    #H_leak = 6.14e4 * (Hlumen - Hstroma)
    # Here, we convert d_protons_in into a "concentration" by dividing by the volumen
    #d_protons_leak = 6.14e4 * (Hlumen*(np.exp(39.5*Dy)) - Hstroma*np.exp(-39.5*Dy))
    #proton leak rate calculated based on P = 2 x 10^-5 cm/s ==> 6.14 e4 per PSII per s
    #39.5 = F/RT, it seems the H_leak has a relatively small impact as claimed by
    #Mordechay SchGnfeld and Hedva Schickler, FEBS letter 1983
    #The permeability of the thylakoid membrane for protons
    dHin = (net_protons_in - v_KEA - v_CLCE)*lumen_protons_per_turnover
    #It looks like earlier code did not calculate v_KEA into H+ concentrations from numbers
    #v_KEA should be numbers of ion/s across KEA. as indicated in dKlumen
    
    # Here we calculate the change in lumen pH by dividing dHin by the buffering capacity
    dpHlumen= -1*dHin / buffering_capacity 

    dHstroma = 0#(net_protons_stroma + v_KEA + v_CLCE)*lumen_protons_per_turnover/10
    #Assuming the volume of stroma is ten times as that of lumen
    dpHstroma = -1*dHstroma / buffering_capacity
    #***************************************************************************************
    #Calculation of Dy considering all ion movements and thylakoid membrane capatitance
    #***************************************************************************************
    delta_charges=charges_from_PSII+PSI_charge_separations + charges_from_b6f \
                    + d_charge_NDH - v_K_channel - d_H_ATP_or_passive - v_VCCN1-3*v_CLCE
    #This net_Klumen does not represent the total charge caused by K+ movement
    #K+ movement only impacts charges from v_K_channel(added variable in this function)            
    #delta_charges= net_protons_in + net_Klumen # - PSII_recombination_v 
    # recall that PSII_recnotesombination is negative electrogenic 
    #note, I now inclluded this term in the calculation of PSII charge separations
    
    dDy=delta_charges*Volts_per_charge
    dpmf= 0.06* dpHlumen + dDy

    #calculate changes to deltaGatp
    #assume that deltaGatp is constant (as suggested by past resarch)...is this changes, 
    #need to consider many metabilic reactions as well.
    #Here we try to consider CBC only
    #DeltaGatp = 30.0 + 2.44* np.log(ATP_pool/ADP_pool/Pi)

    #ddeltaGatp = deltaGatp - DeltaGatp
    #d_ATP_consumed = NADPH_pool*k_CBC*(1-np.exp(-t/900))*1.5
    #if d_ATP_made - d_ATP_consumed < 0:
    #    dATP_pool = 0
    #else:
    dATP_pool= d_ATP_made - d_ATP_consumed
    dADP_pool= - dATP_pool
    #calculate changes in the concentrations of zeaxanthin (Z) and violaxanthin (V)
    #considering VDE_max_turnover_number, pKvde, VDE_Hill, kZE, and lumen pH
    
    dZ, dV = calc_v_VDE(VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pHlumen, V, Z)

    #***************************************************************************************
    #The following calculated changes in NPQ based on the previous and new lumen pH
    #***************************************************************************************

    #calculate the protonation state of PsbS, considering 
    #its pKa and lumen pH
    
    new_PsbS_H = calc_PsbS_Protonation(pKPsbS, pHlumen + dpHlumen)
    new_Z=Z+dZ
    
    #calculate NPQ, based on a simple relationahip between
    #the concentration of Z and the protonation state of PsbS
    #Half contribution from Z but mostly PsbS dependent, half from PsbS alone
    new_NPQ=0.4*max_NPQ*new_PsbS_H*new_Z+0.5*max_NPQ*new_PsbS_H+0.1*max_NPQ*new_Z
    
    #feed this into odeint by calculating the change in NPQ compared to the previous
    #time point
    dNPQ=new_NPQ-NPQ #new_PsbS_H-PsbS_H

    #we re-calculate Phi2 at the start of each iteration of f, so we do not want 
    #odeint to change it
    dPhi2=0 #
    #dADP_pool= 0
    #dATP_pool = 0
    ddeltaGatp = 0
    return [dQA, dQAm, dPQ, dPQH2, dHin, dpHlumen, dDy, dpmf, ddeltaGatp, dKlumen, dKstroma, 
            d_ATP_made, d_PC_ox, d_PC_red, d_P700_ox, d_P700_red, dZ, dV, dNPQ, dsingletO2, dPhi2, dLEF, 
            dFd_ox, dFd_red,  dATP_pool, dADP_pool, dNADPH_pool,dNADP_pool, dCl_lumen, dCl_stroma,dHstroma, dpHstroma]
            

#log_progress displays a tiem bar so the users know how long they have to wait
#for thersults.

def log_progress(sequence, every=None, size=None):

    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = int(size / 200)     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{index} / ?'.format(index=index)
                else:
                    progress.value = index
                    label.value = u'{index} / {size}'.format(
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = str(index or '?')
        
#do_comple_sim does just what it's name says.
#The user sends the initial states (y), the 
#the constants_set_and_trace_times, whcih contains the timing and 
#light intensity (and in the future other parameters),
#and the set of constants (Kx) to use for the simulaitons


def sim(K, initial_states, pulse_times_and_light, max_light_change=1, points_per_segment=1000, **keyword_parameters):
    
    if ('dark_equilibration' in keyword_parameters):
        equibrate_time= keyword_parameters['dark_equilibration']
        use_initial_states=dark_equibration(initial_states, K, equibrate_time)
    else:
        use_initial_states=initial_states
    sub_arrays_time_and_light= optimized_time_split(pulse_times_and_light, 
                                max_light_change, points_per_segment) 

    #first make the constants_set and trace_times for the split segments

    constants_set_and_trace_times=make_variable_light_constants_set_and_trace_times(K, sub_arrays_time_and_light)

    #next, do the simulation
    output=do_complete_sim(use_initial_states, constants_set_and_trace_times, K)
    return(output, use_initial_states)

def sim_ivp(K, initial_states, t_end):    
    output=do_complete_sim(initial_states, t_end, K)
    return output

def do_complete_sim(y00, t_end, Kx):
    
    y00[11]=0 #set ATP_made to zero
    
    # prepare a dictionary with empty arrays to store output
    output={}
    for label in species_labels:
        output[label] = []
    
    #Iterate through the constants_set list, one set of values for each 
    #subtrace, during which the light intensity is held constant    
        
    soln = solve_ivp(f, [0, t_end], y00, args=Kx.as_tuple(), method = 'BDF', \
                     t_eval = np.linspace(0, t_end, 10*t_end+1), max_step = 5)
    
    # Fix the problem with zero time difference between runs.

    time_axis = soln.t
            
        #append a set of computed constants to the output arrays
    for index, label in enumerate( species_labels ):
        output[label] = np.append( output[label], soln.y[index,:] )

    # save the results in case we want to start another simulation that starts where this one left off
    end_state = list(soln.y[:,-1])

    #The following section calculates a number of new parameters from the simulaiton data    
    # Calculate pmf_total from Dy and delta_pH
    Dy = output['Dy']
    pHlumen = output['pHlumen']
    pHstroma = output['pHstroma']
    pmf_total= Dy + ((pHstroma-pHlumen)*.06)
    
    # calculate the Phi2 values based on simulation output parameters

    Phi2_array=[] #contains the calculated Phi2 results
    QA = output['QA']
    NPQ_array = output['NPQ_array']
    for i in range(len(QA)):
        Phi2_array.append(Calc_Phi2(QA[i], NPQ_array[i]))
        
    #calculate tPhiNO and PhiNPQ
    # using the Calc_PhiNO_PhiNPQ function.
    PhiNPQ_array=[]
    PhiNO_array=[]
    for i in range(len(QA)):
        PhiNO, PhiNPQ=Calc_PhiNO_PhiNPQ(Phi2_array[i], QA[i], NPQ_array[i])
        PhiNPQ_array.append(PhiNPQ)
        PhiNO_array.append(PhiNO)
    output['PhiNPQ']=PhiNPQ_array
    output['PhiNO']=PhiNO_array

    #Set up an array to contain the light curve (the PAR values), 
    light_curve=[]    
    for a_t in time_axis:
        light_curve.append(light(a_t, 1200, LIGHT, FREQUENCY, 900, 100))

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
    
    # Some of the values are 
    # duplicates of existing results, with different keys. 
    # Shouldn't be necessary, but I'm leaving this in for now because
    # other function may be expecting these keys
    
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
    
    #calculate the rate of ATP formation, by taking the derivative of the total ATP
    #accumulation
    
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
    
    #output the PsbS protonation state, and the control of electron flow
    #at the cytochrome b6f complex.
    
    output['PsbS_protonated']=[]
    output['b6f_control']=[]
    for pH in output['pHlumen']:
        output['PsbS_protonated'].append(calc_PsbS_Protonation(Kx.pKPsbS, pH))
        output['b6f_control'].append(calc_k_b6f(Kx.max_b6f,1, pH, Kx.pKreg))
        
    #Calculate and store the rate of Fd reduction
    Fd_rate=[0]
    for index in range(1,len(output['Fd_red'])):
        Fd_rate.append((output['Fd_red'][index]-output['Fd_red'][index-1])/(output['time_axis'][index]-output['time_axis'][index-1]))
    output['Fd_rate']=np.array(Fd_rate)
    output['ATP/NADPH']= 2*output['ATP_rate']/(output['Fd_rate'])
    
    ##calculate and store the net fluxes of the counterion K+ into the lumen.
    #K_flux=[0]
    #for i in range(1,len(output['Klumen'])):
    #    K_flux.append((output['Klumen'][i-1]-output['Klumen'][i])/(output['time_axis'][i]-output['time_axis'][i-1]))
    #    
    #output['K_flux']=np.array(K_flux)
    #output['K_flux_normalized']=output['K_flux']/Kx.lumen_protons_per_turnover
    
    #calculate the fluxes of counter-ions and the ratio of LEF_to_NADPH production
    K_flux=[0] #start with zero because we will be calculating the derivative
    for i in range(1,len(output['Klumen'])):
        K_flux.append((output['Klumen'][i-1]-output['Klumen'][i])/(output['time_axis'][i]-output['time_axis'][i-1]))

    output['K_flux']=np.array(K_flux)
    output['K_flux']=output['K_flux']/Kx.lumen_protons_per_turnover
    
    
    # output['LEF_productive']=np.array(output['LEF_productive'])
    # Eliminate nans in the ratio calculations. These occur when flux is zero, when the 
    # ratio is undefined
    for i in range(len(output['ATP_rate'])):
        if np.isnan(output['LEF_to_NADPH'][i]):
            output['LEF_to_NADPH'][i]=output['LEF_to_NADPH'][i-1]
        if np.isnan(output['ATP_rate'][i]):
            output['ATP_rate'][i]=output['ATP_rate'][i-1]
        if np.isnan(output['K_flux'][i]):
            output['K_flux'][i]=output['K_flux'][i-1]

    # calculate the deficit in ATP/NADPH and store it in output['deficit']
    output['deficit']=(output['LEF_to_NADPH']*(3.0/Kx.n)-output['ATP_rate'])  
    output['deficit_int']=integrate.cumtrapz(output['deficit'], output['time_axis'], initial=0)
    output['fract_deficit']=output['deficit_int']/output['LEF_to_NADPH']

    return(output)

# The following coducts a dark equilibration simulation, to allow the system to achieve 
# a steady, dark conditions. 

def dark_equibration(y_initial, Kx, total_duration, **keyword_parameters): 
    #make a sin wave with zero amplitude
    light_frequency=1/total_duration
    points_per_second=10
    max_PAR=0
    dark_time_light_profile=generate_sin_wave(total_duration, max_PAR, light_frequency, points_per_second)
    max_light_change=10
    points_per_segment=100
    optimized_dark_sub_arrays= optimized_time_split(dark_time_light_profile, 
        max_light_change, points_per_segment) 

    #Generate the constants_set and trace_times for the split segments
    constants_set_and_times=make_variable_light_constants_set_and_trace_times(Kx, optimized_dark_sub_arrays)
    
    #next, do the simulation
    output=do_complete_sim(y_initial, constants_set_and_times, Kx)

    #store the final state in  dark_equilibrated_initial_y
    dark_equilibrated_initial_y=output['end_state']

    if ('return_kinetics' in keyword_parameters) and keyword_parameters['return_kinetics']==True:
        return(dark_equilibrated_initial_y, output)
    else:
        return(dark_equilibrated_initial_y)


#*******************************************************************************
#*******************************************************************************
#                   Code related to plotting out results                       *
#*******************************************************************************
#*******************************************************************************


# The plot_interesting_stuff function plots out several graphs
# that show interesting or important results
# It is not meant to output final results

def plot_interesting_stuff(figure_name, output):
    #matplotlib.rcParams.update['figure.figsize'] = [10, 8]
    #light=output['light_curve']
    ltc='red'
    plt.rcParams.update({'font.size': 5})
    time_axis_seconds=output['time_axis']
    max_time=np.max(time_axis_seconds)
    if max_time/(60*60)>1:    
        time_axis=time_axis_seconds/(60*60)
        time_label='Time (h)'
    elif max_time/(60)>1:
        time_axis=time_axis_seconds/(60)
        time_label='Time (min)'
    else:
        time_axis=time_axis_seconds
        time_label='Time (s)'

    fig = plt.figure(num=figure_name, figsize=(5,4), dpi=200)
    ax1 = fig.add_subplot(331)
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=True))
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax1b = ax1.twinx()
    ax1.plot(time_axis, output['pmf'], label='pmf', zorder=3)
    ax1b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax1b.fill_between(time_axis,output['light_curve'],0,color=ltc, alpha=.1, zorder=2)
    ax1.set_xlabel(time_label)
    ax1.set_ylabel('pmf (V)')
    ax1b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax1.set_xlim(0, 1.1*np.max(time_axis))
    ax1b.set_xlim(0, 1.1*np.max(time_axis))
    ax1b.set_ylabel('intensity')
    ax1.yaxis.label.set_color('blue')
    ax1b.yaxis.label.set_color(ltc)
    ax1.spines['left'].set_color('blue')
    ax1b.spines['right'].set_color(ltc)
    ax1b.spines['left'].set_color('blue')
    ax1.tick_params(axis='y', colors='blue')
    ax1b.tick_params(axis='y', colors=ltc)
    ax1.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax1.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax1b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax1b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax2 = fig.add_subplot(332)
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax2b = ax2.twinx()
    ax2.plot(time_axis, output['pHlumen'], label='lumen pH')
    ax2.plot(time_axis, output['pHstroma'],color = 'green')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.set_xlabel(time_label)
    ax2b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax2b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax2b.set_ylabel('intensity')
    ax2.set_ylabel('pH of lumen and stroma')
    ax2.set_xlim(0, 1.1*np.max(time_axis))
    ax2b.set_xlim(0, 1.1*np.max(time_axis))
    ax2b.yaxis.set_major_formatter(FormatStrFormatter('%.f'))
    ax2b.set_ylabel('intensity')
    ax2.yaxis.label.set_color('blue')
    ax2b.yaxis.label.set_color(ltc)
    ax2.spines['left'].set_color('blue')
    ax2b.spines['right'].set_color(ltc)
    ax2b.spines['left'].set_color('blue')
    ax2.tick_params(axis='y', colors='blue')
    ax2b.tick_params(axis='y', colors=ltc)
    ax2.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax2.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax2b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax2b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    
    ax3 = fig.add_subplot(333)
    ax3.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax3b = ax3.twinx()
    ax3.plot(time_axis, output['Dy'], label='Dy')
    ax3b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax3.set_xlabel(time_label)
    ax3.set_ylabel(r'$\Delta\psi$ (V)')
    ax3b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax3b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax3b.set_ylabel('intensity')
    ax3.set_xlim(0, 1.1*np.max(time_axis))
    ax3b.set_xlim(0, 1.1*np.max(time_axis))
    ax3b.set_ylabel('intensity')
    ax3.yaxis.label.set_color('blue')
    ax3b.yaxis.label.set_color(ltc)
    ax3.spines['left'].set_color('blue')
    ax3b.spines['right'].set_color(ltc)
    ax3b.spines['left'].set_color('blue')
    ax3.tick_params(axis='y', colors='blue')
    ax3b.tick_params(axis='y', colors=ltc)
    ax3.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax3.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax3b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax3b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    

    ax4 = fig.add_subplot(334)
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=True)
    ax4.yaxis.set_major_formatter(y_formatter)
    ax4.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax4b = ax4.twinx()
    ax4b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax4.plot(time_axis, output['Klumen'], label='K+ lumen')
    ax4.set_xlabel(time_label)
    ax4.set_ylabel(r'$K^{+} in lumen$')
    ax4b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax4b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax4.set_xlim(0, 1.1*np.max(time_axis))
    ax4b.set_xlim(0, 1.1*np.max(time_axis))
    ax4b.set_ylabel('intensity')
    ax4.yaxis.label.set_color('blue')
    ax4b.yaxis.label.set_color(ltc)
    ax4.spines['left'].set_color('blue')
    ax4b.spines['right'].set_color(ltc)
    ax4b.spines['left'].set_color('blue')
    ax4.tick_params(axis='y', colors='blue')
    ax4b.tick_params(axis='y', colors=ltc)
    ax4.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax4.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax4b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax4b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 

    ax5 = fig.add_subplot(335)
    ax5.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax5b = ax5.twinx()
    ax5b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax5.plot(time_axis, 1-output['QAm'], color='green', label='qL')
    ax5.plot(time_axis, output['Phi2'], color = 'blue',label ='Phi2')
    ax5b.plot(time_axis, output['NADPH_pool'], color='red', label='NADPH_pool')
    ax5.plot(time_axis, output['P700_red'], color = 'm', label = 'P700_red')
    ax5.set_xlabel(time_label)
    ax5.set_ylabel('qL(green) and Phi2')
    ax5b.set_ylabel(r'NADPH_pool')
    ax5.set_xlim(0, 1.1*np.max(time_axis))
    ax5b.set_xlim(0, 1.1*np.max(time_axis))
    ax5.tick_params(axis='y', colors='blue')
    ax5b.tick_params(axis='y', colors='red')
    ax5.yaxis.label.set_color('blue')
    ax5b.yaxis.label.set_color('red')
    ax5.spines['left'].set_color('blue')
    ax5b.spines['right'].set_color('red')
    ax5b.spines['left'].set_color('blue')
    ax5.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax5.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax5b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax5b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 

    ax6 = fig.add_subplot(336)
    ax6.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax6b = ax6.twinx()
    ax6.plot(time_axis, output['QAm'], color='blue', label='QA-')
    ax6b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax6b.plot(time_axis, output['PQH2'], color='green', label='P700_ox')
    ax6.set_xlabel(time_label)
    ax6.set_ylabel(r'$Q_A^{-}$')
    ax6b.set_ylabel(r'$PQH_2$')
    ax6.set_xlim(0, 1.1*np.max(time_axis))
    ax6b.set_xlim(0, 1.1*np.max(time_axis))
    ax6.tick_params(axis='y', colors='blue')
    ax6b.tick_params(axis='y', colors='green')
    ax6.yaxis.label.set_color('blue')
    ax6b.yaxis.label.set_color('green')
    ax6.spines['left'].set_color('blue')
    ax6b.spines['right'].set_color('green')
    ax6b.spines['left'].set_color('blue')
    ax6.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax6.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax6b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax6b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 

    ax7 = fig.add_subplot(337)
    ax7.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax7b = ax7.twinx()
    ax7b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax7.plot(time_axis, output['Z'], label='Z')
    ax7.plot(time_axis, output['V'], label='V')
    ax7.set_xlabel(time_label)
    ax7.set_ylabel('Z and V')
    ax7b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax7b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax7b.set_ylabel('intensity')
    ax7.set_xlim(0, 1.1*np.max(time_axis))
    ax7b.set_xlim(0, 1.1*np.max(time_axis))
    ax7b.set_ylabel('')
    ax7.yaxis.label.set_color('blue')
    ax7b.yaxis.label.set_color(ltc)
    ax7.spines['left'].set_color('blue')
    ax7b.spines['right'].set_color(ltc)
    ax7b.spines['left'].set_color('blue')
    ax7.tick_params(axis='y', colors='blue')
    ax7b.tick_params(axis='y', colors=ltc)
    ax7.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax7.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax7b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax7b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax7c = ax7.twinx()

    #ax7.plot(time_axis, output['V'], label='V')
    ax7c.spines['right'].set_color('orange')
    ax7c.plot(time_axis, output['PsbS_protonated'], label='PsbSH+', color='orange')
    ax7c.set_ylabel('PsbSH+')
    ax7c.yaxis.label.set_color('orange')
            
    ax8 = fig.add_subplot(338)
    ax8.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax8b = ax8.twinx()
    ax8b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax8.plot(time_axis, output['NPQ'], label='qE')
    ax8.set_xlabel(time_label)
    ax8.set_ylabel('NPQ (qE)')
    ax8b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax8b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax8.set_xlim(0, 1.1*np.max(time_axis))
    ax8b.set_xlim(0, 1.1*np.max(time_axis))
    ax8b.set_ylabel('intensity')
    ax8b.set_ylabel('intensity')
    ax8.yaxis.label.set_color('blue')
    ax8b.yaxis.label.set_color(ltc)
    ax8.spines['left'].set_color('blue')
    ax8b.spines['right'].set_color(ltc)
    ax8b.spines['left'].set_color('blue')
    ax8.tick_params(axis='y', colors='blue')
    ax8b.tick_params(axis='y', colors=ltc)
    ax8.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax8.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax8b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax8b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 

    ax9 = fig.add_subplot(339)
    ax9.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax9b = ax9.twinx()
    ax9b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax9.plot(time_axis, output['Cl_lumen'], color='blue', label='Cl_lumen')
    ax9.plot(time_axis, output['Cl_stroma'],color = 'red', label ='Cl_stroma')
    ax9.set_xlabel(time_label)
    ax9.set_ylabel(r'$Cl^{-} lumen(blue), stroma(red)$')
    ax9b.fill_between(time_axis,output['light_curve'],0,color='red', alpha=.1)
    ax9b.set_ylim(0, 1.1*np.max(output['light_curve']))
    ax9.set_xlim(0, 1.1*np.max(time_axis))
    ax9b.set_xlim(0, 1.1*np.max(time_axis))
    ax9.yaxis.label.set_color('blue')
    ax9b.yaxis.label.set_color(ltc)
    ax9.spines['left'].set_color('blue')
    ax9b.spines['right'].set_color(ltc)
    ax9b.spines['left'].set_color('blue')
    ax9.tick_params(axis='y', colors='blue')
    ax9b.tick_params(axis='y', colors=ltc)
    ax9.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax9.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 
    ax9b.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
    ax9b.locator_params(axis = 'y', nbins = 4)# (or axis = 'y') 

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.show()



#                    
#find the global max and min for a list of arrays

def global_min_max(list_of_arrays):
    local_max=[]
    local_min=[]
    for local_array in list_of_arrays:
        local_min.append(np.min(np.array(local_array)))
        local_max.append(np.max(np.array(local_array)))
        #print(local_max)
    global_min=np.min(local_min)
    global_max=np.max(local_max)
    return (global_min,global_max)
                    

def get_axis_limits(ax, scale=.9):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale

# plot_gen is a generalized routine for plotting out the kinds of graphs I use
#for the simulation data
# fig = the figure object to which to plot
# sub_plot_number = the subplot number to which to add the plot data
# plot_list is the list of x,y and otehr parameters
# plot_every_nth_point will tell the plotting routine to ony plot a certain
# number of points.
# More details in the code:

def plot_gen(fig, sub_plot_number, plot_list, plot_every_nth_point, **keyword_parameters):
        
    #make three axes, two for data, one for light curve(if needed)
    #all have same x-axis

    any_left=False
    any_right=False
    any_light=False

    for plots in plot_list:
        #print(plots.what_to_plot[1])
        if plots.axis == 'left':
            any_left=True
        if plots.axis == 'right':
            any_right=True
        if plots.axis == 'light':
            any_light=True

    all_axes=[]
    ax1a=fig.add_subplot(sub_plot_number[0], sub_plot_number[1], sub_plot_number[2]) #I assume we have a left axis graph
    
    ax1a.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    all_axes.append(ax1a)
    if any_right: #if we have any right axis graphs, add axis
        ax1b= ax1a.twinx()
        all_axes.append(ax1b)
        ax1b.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
    if any_light: #if we have any light axis graphs, add axis
        #print('found light curve')
        ax1c = ax1a.twinx()
        all_axes.append(ax1c)
        ax1c.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    
    for plots in plot_list: #iterate through all the things we want to plot
        output=plots.output

        #the following makes invisible, the inside facing axes. Only the left-most left 
        #axis and the rigth-most right axis lables will be shown.          
        #say we have a 3 conditions and 2 phenomena. 
        # We want to make the left axes on the leftmost panels to be visible. The left most panels are:
        # 1, 4 
        #which is when sub_plot_number[2]+ sub_plot_number[1]-1 is divisible by the number of conditons, e.g.
        # 1+3-1 = 3 or 4+3-1 =6...
        # I test for this using the modulus function %
        
        
        
        if (sub_plot_number[2]+ sub_plot_number[1]-1)%sub_plot_number[1]==0 and plots.axis=='left': 
            plots.yaxis_visible=True
            
            #next we check for the rightmost panels, noting that in this case, the panel numbers are:
            #3 and #6, both of which are divisible by 3 
            
        elif int(sub_plot_number[2])%int(sub_plot_number[1])==0 and plots.axis=='right':
                plots.yaxis_visible=True
        else:
            # if neither of these conditions are true, we make the axes invisible
            plots.yaxis_visible=False
            plots.show_legend=False  #we only want one copy of the legend. If plots.show_gend is True,
                                     #it will only show for the left most or right-most panels 

        if plots.axis=='left': #use axis ax1a
            ax1=ax1a
            ax1.yaxis.label.set_color(plots.axis_color)
            ax1.spines['left'].set_color(plots.axis_color)
            ax1.tick_params(axis='y', colors=plots.axis_color)
            ax1.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
            ax1.locator_params(axis = 'y', nbins = 6)# (or axis = 'y') 
            plot_font_size=plots.plot_font_size
            #print(plot_font_size)
            ax1.set_xlabel(plots.x_axis_label, fontsize=plot_font_size*1.1)
            ax1.set_ylabel(plots.y_axis_label, fontsize=plot_font_size*1.1, labelpad=2)
            #ax1.set_xlim(0, 1.2*np.max(x_values))
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), fontsize=plot_font_size*1.2)

        elif plots.axis=='right':
            ax1=ax1b
            ax1.yaxis.label.set_color(plots.axis_color)
            ax1.spines['right'].set_color(plots.axis_color)
            ax1.tick_params(axis='y', colors=plots.axis_color)
            ax1.locator_params(axis = 'x', nbins = 4)# (or axis = 'y') 
            ax1.locator_params(axis = 'y', nbins = 6)# (or axis = 'y') 
            plot_font_size=plots.plot_font_size
            ax1.set_xlabel(plots.x_axis_label, fontsize=plot_font_size)
            ax1.set_ylabel(plots.y_axis_label, size=plot_font_size*1.2, labelpad=2)
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), fontsize=plot_font_size*1.2)
            
        elif plots.axis=='light':
            ax1=ax1c
            ax1.spines['right'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.axes.get_yaxis().set_visible(False)
        
        #detect the bottom row, all others set x-axis invisible
        
        if (sub_plot_number[1] * (sub_plot_number[0]-1)) > sub_plot_number[2]-1:
            ax1.set_xlabel('')
            plt.setp(ax1.get_xticklabels(), visible=False)
        
        x_values=output[plots.what_to_plot[0]][::plot_every_nth_point]

        if plots.subtract_baseline==True:
            y_values= output[plots.what_to_plot[1]]-output[plots.what_to_plot[1]][0]
            y_values=y_values[::plot_every_nth_point]
        else:
            y_values= output[plots.what_to_plot[1]][::plot_every_nth_point]

        ax1.ticklabel_format(useOffset=False)

        if plots.zero_line==True:
            zero_line=[np.min(x_values),np.max(x_values)]
            
            ax1.plot(zero_line,[0,0], linestyle='--', color='grey')

        if plots.axis=='light':
            ax1.fill_between(x_values, y_values,0,color=plots.marker_color, alpha=.1)
        else:
            if plots.linestyle=='solid':                   
                ax1.plot(x_values, y_values, color=plots.marker_color, label=plots.data_label, 
                        lw=plots.linewidth, zorder=3, linestyle=plots.linestyle)
            elif plots.linestyle=='dashed':
                ax1.plot(x_values, y_values, color=plots.marker_color, label=plots.data_label, 
                        lw=plots.linewidth, zorder=3, linestyle=plots.linestyle, dashes=(3, 1))
            else:
                ax1.plot(x_values, y_values, color=plots.marker_color, label=plots.data_label, 
                        lw=plots.linewidth, zorder=3, linestyle=plots.linestyle, dashes=(1, 2))

        if plots.set_maxmin_y==True:
            ax1.set_ylim(plots.maxmin_y[0], plots.maxmin_y[1])  
        else:
            
            ypad=0.1*(np.max(y_values)-np.min(y_values))
            
            ax1.set_ylim(np.min(y_values)-ypad, np.max(y_values)+ypad)  
            
        if plots.set_maxmin_x==True:
            ax1.set_xlim(plots.maxmin_x[0], plots.maxmin_x[1])  
        else:
            ax1.set_xlim(np.min(x_values), np.max(x_values))  
        if plots.axis=='light':
            ax1.set_ylim(0, np.max(y_values))  

        if plots.show_legend==True:
            ax1.legend(loc='upper center', bbox_to_anchor=(0.75, 0.99), fancybox=False, 
                       shadow=False, frameon=False, ncol=1,
                      fontsize=plot_font_size)
        if plots.yaxis_visible==False:
                ax1.axes.get_yaxis().set_visible(False)
        else:
                ax1.axes.get_yaxis().set_visible(True)
        sub_plot_annotation=''
        if ('annotation' in keyword_parameters):
            sub_plot_annotation=keyword_parameters['annotation']
        # place a text box in upper left in axes coords
        props = dict(boxstyle='circle', facecolor='white', alpha=1)
        ax1.text(0.8, .2, sub_plot_annotation, transform=ax1.transAxes, fontsize=plot_font_size,
                    verticalalignment='top', bbox=props)
        ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2), fontsize=plot_font_size)
        ax1.tick_params(labelsize=plot_font_size)
    return(all_axes)





def plot_pmf_params(output, use_x_axis, x_axis_label, global_min, global_max):
    sub_base=False
    
    all_min=np.min([global_min['Dy'], global_min['delta_pH_V'],global_min['pmf'] ])
    all_max=np.max([global_max['Dy'], global_max['delta_pH_V'],global_max['pmf'] ])
    
    #set up the left axis of the plot for membrane potential
    what_to_plot=[use_x_axis, 'Dy']
    a1=sim_plot()
    a1.output=output
    a1.what_to_plot=what_to_plot
    a1.data_label=r'$\Delta\psi$'
    a1.axis_color='black'
    a1.marker_color='blue'
    a1.linestyle='solid'
    a1.y_axis_label='potential (V)'
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label
    else:
        a1.x_axis_label='time (s)'
    a1.axis='left'
    a1.linewidth=1
    a1.subtract_baseline=sub_base
    a1.plot_font_size=7
    a1.zero_line=True
    a1.set_maxmin_y=True
    a1.set_maxmin_x=False
    a1.maxmin_x=[0,1000]
    
    a1.maxmin_y=[all_min, all_max]

    a1.yaxis_visible=True
    a1.show_legend=True

    #add to the left axis, a plot for delta_pH 

    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'delta_pH_V']
    #plot delta pH
    a2.data_label=r'$\Delta$pH'
    a2.axis_color='black'
    a2.marker_color='red'
    a2.what_to_plot=what_to_plot
    a2.linestyle='dashed'
    a2.subtract_baseline=sub_base
    a2.maxmin_y=[all_min, all_max]


    a3=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'pmf']
                            
    # add to the left axis, a plot of pmf in V
    a3.what_to_plot=what_to_plot
    a3.data_label='pmf'
    a3.axis_color='black'
    a3.marker_color='green'
    a3.linestyle='solid'
    a3.subtract_baseline=sub_base

    #a3.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a3.maxmin_y=[all_min, all_max]



    #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=False
    #a4.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a4.maxmin_y=[all_min, all_max]
    a4.set_maxmin_x=False
    a4.yaxis_visible=False
    a4.maxmin_x=[0,1000]
    
    
    return [a1,a2,a3,a4]

def plot_pmf_params_offset(output, use_x_axis, x_axis_label, global_min, global_max):
    sub_base=False
    
    all_min=np.min([global_min['Dy_offset'], global_min['delta_pH_V_offset'],global_min['pmf_offset'] ])
    all_max=np.max([global_max['Dy_offset'], global_max['delta_pH_V_offset'],global_max['pmf_offset'] ])
    
    #set up the left axis of the plot for membrane potential
    what_to_plot=[use_x_axis, 'Dy_offset']
    a1=sim_plot()
    a1.output=output
    a1.what_to_plot=what_to_plot
    a1.data_label=r'$\Delta\psi$'
    a1.axis_color='black'
    a1.marker_color='blue'
    a1.linestyle='solid'
    a1.y_axis_label=r'$\Delta$ V'
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label
    else:
        a1.x_axis_label='time (s)'
    a1.axis='left'
    a1.linewidth=1
    a1.subtract_baseline=sub_base
    a1.plot_font_size=7
    a1.zero_line=True
    a1.set_maxmin_y=True
    a1.set_maxmin_x=False
    a1.maxmin_x=[0,1000]
    a1.maxmin_y=[all_min, all_max]
    a1.yaxis_visible=True
    a1.show_legend=True

    #pmf_parameters_plot.append(a1)

    #add to the left axis, a plot for delta_pH 

    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'delta_pH_V_offset']
    #plot delta pH
    a2.data_label=r'$\Delta$pH'
    a2.axis_color='black'
    a2.marker_color='red'
    a2.what_to_plot=what_to_plot
    a2.linestyle='dashed'
    a2.subtract_baseline=sub_base
    a2.maxmin_y=[all_min, all_max]


    a3=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'pmf_offset']
                            
    # add to the left axis, a plot of pmf in V
    a3.what_to_plot=what_to_plot
    a3.data_label='pmf'
    a3.axis_color='black'
    a3.marker_color='green'
    a3.linewidth=1.5
    a3.linestyle='dashed'
    a3.subtract_baseline=sub_base
    a3.maxmin_y=[all_min, all_max]


    #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=False
    a4.maxmin_y=[all_min, all_max]
    a4.set_maxmin_x=False
    a4.yaxis_visible=False
    a4.maxmin_x=[0,1000]
    return [a1,a2,a3,a4]
    

def plot_K_and_parsing(output, use_x_axis, x_axis_label,global_min, global_max):

    #set up the left axis of the plot for NPQ
    
    a1=sim_plot() #define an instance of the plot 
    what_to_plot=[use_x_axis, 'NPQ'] #indicate thwat to plot. the variable 'use_this_axis' is passed to the function
    
    a1.y_axis_label=r'q$_{E}$' # set the y-axis label
    a1.data_label='qE'
    a1.output=output
    
    a1.what_to_plot=what_to_plot
    a1.axis_color='green'
    a1.marker_color='green'
    a1.linestyle='dashed'
    a1.y_axis_label='NPQ (qE)'
    
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label #if there is something in x_axis_label then use it.
    else:
        a1.x_axis_label='time (s)'
    a1.axis='left'
    a1.linewidth=1
    a1.subtract_baseline=False
    a1.plot_font_size=7
    a1.zero_line=False
    a1.set_maxmin_y=True
    a1.set_maxmin_x=False
    a1.maxmin_x=[0,1000]
    a1.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a1.yaxis_visible=True
    a1.show_legend=False

    #set up the right axis of the plot for [K+]

    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'Klumen']
    a2.data_label=r'lumen $K^{+}$ (M)' #'[K+] lumen (M)'
    a2.axis_color='black'
    a2.marker_color='black'
    a2.what_to_plot=what_to_plot
    a2.linestyle='solid'
    a2.axis='right'
    a2.zero_line=False
    a2.y_axis_label=r'lumen $K^{+}$ (M)' #'[K+] lumen (M)'
    a2.show_legend=False
    a2.set_maxmin_y=True
    #a2.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1][1]]]
    a2.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]

    #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=True
    a4.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a4.set_maxmin_x=False
    a4.yaxis_visible=False
    a4.maxmin_x=[0,1000]
    a4.show_legend=False
    return [a1,a2,a4]
    

        
def b6f_and_balance(output, use_x_axis, x_axis_label, global_min, global_max):

    #set up the left axis of the plot for NPQ

    a1=sim_plot()
    a1.output=output

    what_to_plot=[use_x_axis, 'b6f_control']
    a1.data_label='rate constant for b6f' #'[K+] lumen (M)'
    a1.axis_color='blue'
    a1.marker_color='blue'
    a1.what_to_plot=what_to_plot
    a1.linestyle='dashed'
    a1.axis='left'
    a1.zero_line=False
    a1.y_axis_label= r'$b_{6}f$ rate constant $(s^{-1})$'  # r'k_{bf}' #'[K+] lumen (M)'  r'$^{1}O_{2}$ $(s^{-1})$ (cumulative)'
    a1.show_legend=False
    a1.set_maxmin_y=True
    a1.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a1.set_maxmin_x=False
    a1.yaxis_visible=False
    a1.maxmin_x=[0,1000]
    a1.show_legend=False
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label
        #print('over riding x_axis label')
    else:
        a1.x_axis_label='time (s)'


    a1.yaxis_visible=True
    a1.show_legend=False
    
    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'pHlumen']
    a2.data_label=r'lumen pH'
    a2.axis_color='red'
    a2.marker_color='red'
    a2.what_to_plot=what_to_plot
    a2.linestyle='solid'
    a2.axis='right'
    a2.zero_line=False
    a2.y_axis_label=r'lumen pH'
    a2.show_legend=False
    a2.set_maxmin_y=True
    a2.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    
        #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=True
    a4.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a4.set_maxmin_x=False
    a4.yaxis_visible=False
    a4.maxmin_x=[0,1000]
    a4.show_legend=False
    

    return [a1, a2, a4]
    
    
def plot_QAm_and_singletO2(output, use_x_axis, x_axis_label,global_min, global_max):

    #set up the left axis of the plot for NPQ
    a1=sim_plot()
    what_to_plot=[use_x_axis, 'QAm']
    a1.data_label=r'$Q_A^{-}$'
    a1.output=output
    a1.what_to_plot=what_to_plot
    a1.axis_color='green'
    a1.marker_color='green'
    a1.linestyle='dashed'
    a1.y_axis_label=r'$Q_A^{-}$'
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label
    else:
        a1.x_axis_label='time (s)'
    a1.axis='left'
    a1.linewidth=1
    a1.subtract_baseline=False
    a1.plot_font_size=7
    a1.zero_line=False
    a1.set_maxmin_y=True
    a1.set_maxmin_x=False
    a1.maxmin_x=[0,1000]
    a1.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a1.yaxis_visible=True
    a1.show_legend=False

    #set up the right axis of the plot for 1O2

    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'singletO2_rate']
    a2.data_label=r'$^{1}O_{2}$ $(s^{-1})$'
    a2.axis_color='red'
    a2.marker_color='red'
    a2.what_to_plot=what_to_plot
    a2.linestyle='solid'
    a2.axis='right'
    a2.zero_line=False
    a2.y_axis_label=r'$^{1}O_{2} (s^{-1})$'
    a2.show_legend=False
    a2.set_maxmin_y=True
    a2.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]

    #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=True
    a4.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a4.set_maxmin_x=False
    a4.yaxis_visible=False
    a4.maxmin_x=[0,1000]
    a4.show_legend=False
    return [a1,a2,a4]
    
def plot_cum_LEF_singetO2(output, use_x_axis, x_axis_label,global_min, global_max):
    #set up the left axis of the plot for commulative LEF
    a1=sim_plot()
    what_to_plot=[use_x_axis, 'LEF_cumulative']
    a1.data_label='LEF cumulative'
    a1.output=output
    
    a1.what_to_plot=what_to_plot
    a1.axis_color='green'
    a1.marker_color='green'
    a1.linestyle='dashed'
    a1.y_axis_label='LEF cumulative'
    if x_axis_label != '': 
        a1.x_axis_label=x_axis_label
    else:
        a1.x_axis_label='time (s)'
    a1.axis='left'
    a1.linewidth=1
    a1.subtract_baseline=False
    a1.plot_font_size=7
    a1.zero_line=False
    a1.set_maxmin_y=True
    a1.set_maxmin_x=False
    a1.maxmin_x=[0,1000]
    #a1.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a1.maxmin_y=[0,global_max[what_to_plot[1]]]
    a1.yaxis_visible=True
    a1.show_legend=False

    #set up the right axis of the plot for [K+]
    a2=copy.deepcopy(a1) #sim_plot()
    what_to_plot=[use_x_axis, 'singletO2_array']
    a2.data_label=r'$^{1}O_{2}$ $(s^{-1})$ (cumulative)'
    a2.axis_color='red'
    a2.marker_color='red'
    a2.what_to_plot=what_to_plot
    a2.linestyle='solid'
    a2.axis='right'
    a2.zero_line=False
    a2.y_axis_label=r'$^{1}O_{2}$ (cumulative)'
    a2.show_legend=False
    a2.set_maxmin_y=True
    a2.maxmin_y=[0,global_max[what_to_plot[1]]]
    a2.set_maxmin_x=True
    a2.maxmin_x=[global_min[what_to_plot[0]],global_max[what_to_plot[0]]]

    #set up the right axis of the plot for the llight curve as a semi-transparent filled area, do not discply the axis
    a4=copy.deepcopy(a1)
    what_to_plot=[use_x_axis, 'light_curve']
    a4.axis='light'
    a4.what_to_plot=what_to_plot
    a4.axis_color='red'
    a4.marker_color='red'
    a4.subtract_baseline=False
    a4.set_maxmin_y=True
    a4.maxmin_y=[global_min[what_to_plot[1]],global_max[what_to_plot[1]]]
    a4.set_maxmin_x=False
    a4.set_maxmin_x=True
    a4.maxmin_x=[global_min[what_to_plot[0]],global_max[what_to_plot[0]]]
    a4.show_legend=False

    return [a1,a2,a4]
    
def best_time_scale(output):
    #determine the best time axis to use
    max_time=np.max(output['time_axis'])
    if max_time>3599:    
        use_time_axis='time_axis_h'
        time_label='Time (h)'
    elif max_time>60:
        use_time_axis='time_axis_min'
        time_label='Time (min)'
    else:
        use_time_axis='time_axis'
        time_label='Time (s)'
    return(use_time_axis, time_label)
    
def find_global_max_min(output_list, conditions_to_plot, pad):
    global_min={}
    global_max={}
    rep_output=list(output_list.keys()) #pick the first item on output_list as a representative 

    for k in output_list[rep_output[0]]: #iterate through all the data arrays in output
        gmin_array=[]
        gmax_array=[]
        for condition_name in conditions_to_plot:
            #print(condition_name)
            output=output_list[condition_name]
            gmin_array.append(np.min(output[k]))
            gmax_array.append(np.max(output[k]))
        gmin=np.min(gmin_array)
        gmax=np.max(gmax_array)
        global_min[k]=gmin-(pad*(gmax-gmin))
        global_max[k]=gmax+(pad*(gmax-gmin))
    return(global_min, global_max)
    
#generate a dictionary of plotting functions so they can be more easily called in loops

plot_results={}
plot_results['pmf_params']=plot_pmf_params
plot_results['pmf_params_offset']=plot_pmf_params_offset

plot_results['K_and_parsing']=plot_K_and_parsing
plot_results['plot_QAm_and_singletO2']=plot_QAm_and_singletO2
plot_results['plot_cum_LEF_singetO2']=plot_cum_LEF_singetO2
plot_results['b6f_and_balance'] = b6f_and_balance

def plot_block(output_list, fig, conditions_to_plot, where, phenomena_sets, plot_every_nth_point):
    
    subplot_col_labels=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    subplot_row_labels=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    
    #determine the number of colums
    num_cols=len(where)
    num_rows=len(conditions_to_plot)
    number_phenomena_sets=len(phenomena_sets)
    
    global_min, global_max=find_global_max_min(output_list, conditions_to_plot, .1)

    for j, phenomena in enumerate(phenomena_sets):
        for i, condition_name in enumerate(conditions_to_plot):
            output=output_list[condition_name]

            #determine the best time axis (s, min, hours) to use: 
            use_time_axis, time_label=best_time_scale(output)

            #determine the sub_plot_number from number_phenomena_sets, num_rows, j, and where[i]]
            sub_plot_number=[number_phenomena_sets,num_rows,(j*num_rows)+where[i]]
            if j+1==num_rows:
                #print('bottom')
                time_label=''
                
            plot_list=plot_results[phenomena](output, use_time_axis, time_label, global_min, global_max)

            #subplot_annotation_number=j*num_rows+i
            an_text=str(subplot_col_labels[i]) + '.' + subplot_row_labels[j]
            
            plot_gen(fig, sub_plot_number, plot_list, plot_every_nth_point, 
            subplot_label=an_text, annotation=an_text) #subplot_lables[subplot_annotation_number])
            
#shrink will decrease the size of very large data sets by returning every nth point sets

def shrink(output, start_point, take_every_nth_point):
    shrunk={}
    for key in list(output.keys()):
        shrunk[key]=output[key][::take_every_nth_point]
    return(shrunk)


    
class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

class sim_plot(FrozenClass):
    def __init__(self):
        self.linewidth=1
        self.output={}
        self.data_label=r'$\Delta\psi$'
        self.axis='left'
        self.maxmin=[] #[-1,1]
        self.what_to_plot=['time_axis', 'Dy']
        self.axis_color='blue'
        self.marker_color='blue'
        self.subtract_baseline=False
        self.plot_font_size=7
        self.linestyle='solid'
        self.y_axis_label='y axis'
        self.x_axis_label='x axis'
        self.zero_line=True
        self.set_maxmin_y=True
        self.set_maxmin_x=False
        self.maxmin_x=[0,1000]
        self.maxmin_y=[-.05,.27]
        self.yaxis_visible=True
        self.show_legend=True
        self._freeze() # no new attributes after this point.


#Set up STANDARD initial conditions, most of these values were taken from Cruz et al., 2005
class standard_constants(object):
    #***************************************************************************************
    #paramweters for V-->Z and Z-->V reactions 
    #***************************************************************************************
    VDE_max_turnover_number=0.08#changed from 1
    pKvde= 5.65#changed from 6.0
    #"reviewed in: Peter Jahns a,⁎, Dariusz Latowski b,c, Kazimierz StrzalkaMechanism and
    #regulation of the violaxanthin cycle: The role of antenna proteins and
    #membrane lipids. Biochimica et Biophysica Acta 1787 (2009) 3–14" Erhard E. Pfündel*2
    #and Richard A. Dille, The pH Dependence of Violaxanthin Deepoxidation in lsolated
    #Pea Chloroplasts. Plant Physiol. (1993) 101: 65-71

    VDE_Hill=4 
    kZE=0.004 #changed from 0.01

    #***************************************************************************************
    #paramweters for PsbS protonation 
    #***************************************************************************************

    pKPsbS=6.2  #pKa for protonation of PsbS, 6.0-6.5, assuming Hill coefficient=1
    max_NPQ=3  #this max_NPQ is how much PsbS can result in, Zeaxanthin can play half of it
    # but PsbS dependent, PsbS can independent play half of it

    pHstroma_initial=7.8
    #***************************************************************************************
    #paramweters for ATP SYNTHASE
    #***************************************************************************************
    ATP_synthase_content= 0.367 #0.5 for Photosynth Res 2013 117:1-30
    # or 0.367 ##Journal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014
    #doi:10.1093/jxb/eru090 Advance Access publication 12 March, 2014
    #Molecular Architecture of the Thylakoid Membrane: Lipid Diffusion Space for
    #Plastoquinone, H. Kirchhoff,*,| U. Mukherjee,| and H.-J. Galla§ Biochemistry 2002, 41, 4872-4882
    #Heinrich Strotmann and Susanne Bickel-Sandkotter, STRUCTURE, FUNCTION, AND REGULATION OF 
    #CHLOROPLAST ATPase. Ann. Rev. Plant Physiol. 1984. 35:97-120

    ATP_synthase_max_turnover=200.0#per PSII
    #or one can use in vitro data of maxrate /ATPase is 400 then
    #ATP_synthase_max_turnover=400*ATP_synthase_content
    #this would give a max_turnover from 146.8 to 200 per PSII
    
    #However, one may also consider that there is a maximal (saturating turover rate 
    #(saturation point), as shown by Junesch and Grabber (1991)
    #http://dx.doi.org/10.1016/0014-5793(91)81447-G
    #another ref: Ulrike Junesch and Peter Graber, 1987 BBA 275-288
    #Influence of the redox state an dthe activation of the chloroplast ATP synthase...


    #***************************************************************************************
    #Membrane capacitance
    #***************************************************************************************
    Thylakoid_membrane_capacitance = 0.6
    Volts_per_charge=0.047 #thylakoid membrane capacitance = 0.6 uF/cm2

    
    #print('the DGATP is set to: ' + str(DeltaGatp_KJ_per_mol) + ' kJ/mol, which is: ' + str(DeltaGatp_initial) + ' volts')
    #ATP=1/(10**((32-DeltaGatp_initial)/5.7)+1)

    #***************************************************************************************
    # The value of n, calculated from c subunits
    # Here is where we set up the value of n, the ratio of protons/ATP, which I assume is (# of c subunits)/3
    #***************************************************************************************
    
    
    c_subunits_per_ATP_synthase=14
    n=c_subunits_per_ATP_synthase/3 


    #***************************************************************************************
    # Counter ion exchange reactions
    #***************************************************************************************
    #permeability of the thylakoid to K+ 
    perm_K= 150 #if calculated from 3.6*10^-8cm/s, 510 nm2/PSII, then 111/s
    #Hui Lyu and Dusan Lazar Journal of Theoretical Biology 413 (2017) 11-23
    #one can play with this parameter to see how it affects the kinectics.
    
    #VCCN1 conductance, or Cl- permeability through VCCN1 unit: Cl-/M/V/s
    #use k_VCCN1 = 180 if Cl_flux_relative is not used.
    #if Cl_flux_relative function is used, this is rate constant at 0.1 V
    #unit: Cl-/M/s, normalized to per PSII.
    k_VCCN1 = 12
    #k_CLCE is the kinetic constant of CLCE2.
    k_CLCE = 800000
    #***************************************************************************************
    # b6f reactions
    #***************************************************************************************
    b6f_content=0.433#Journal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014
    #doi:10.1093/jxb/eru090 Advance Access publication 12 March, 2014
    max_b6f=300.0
    pKreg=6.2
    #revew: Tikhonov AV 2014 The cytochrome b6f complex at the crossroad of photosynthetic electron transport pathways. Plant Physiol Biochem 81, 163-183
    Em7_PC=0.37
    Em7_PQH2 = 0.11
    
    Em_Fd = -0.42
    k_NDH = 1000.0 # unit s^-1

    #***************************************************************************************
    # Lumen proton bufering reactions
    #***************************************************************************************
    lumen_protons_per_turnover= 0.000587 #the change in molarity with one H+ transferred to 
    #lumen per PSII
    buffering_capacity=0.014

    PSI_antenna_size=0.5 #setting this to the same valus as PSII_antenna_size will imply that 
                        #equal numbers of photons hit the two photosystems

    k_PC_to_P700=5000 #rate constant for oxidation of PC by P700+

    #***************************************************************************************
    # Proton exchange through the KEA3 system
    #***************************************************************************************

    k_KEA=2500000#this would have significant impact, consider low [H+] in lumen

    #***************************************************************************************
    #parameters for PSII reactions 
    #***************************************************************************************
    max_PSII=1     
    PSII_antenna_size=0.5

    kQA=1000  #the rate constant for oxidation of QA- by PQ


    #***************************************************************************************
    #parameters for recombination and singlet O2 production 
    #***************************************************************************************
    k_recomb=0.33
    triplet_yield=.45 #45% of recomnbinations lead to P680 triplets
    triplet_to_singletO2_yield=1 #100% of 3P680 give rise to 1O2


    #***************************************************************************************
    #Light intensity in terms of hits per second to PSII associated antenna 
    #***************************************************************************************

    light_per_L=0

    k_Fd_to_NADP=1000 #rate constant for transfer of electrons from Fd to NADP
    
    k_CBC=60 #max rate constant for the Calvin-Benson cycle in terms of NADPH consumption
    #this number is adjusted based on light intensity.
    

#***************************************************************************************
#***************************************************************************************
# Initial concentrations and parameters
#***************************************************************************************
#***************************************************************************************


class standard_initial_states(object):
    V_initial=1.0
    Z_initial=0.0
    #start with no ATP made
    ATP_made_initial=0
    DeltaGatp_initial =  30.0 + 2.44 * np.log(1/0.0015)#KJ/mol
    Klumen_initial=0.100
    Kstroma_initial=0.100
    
    Cl_lumen_initial = 0.04
    Cl_stroma_initial = 0.04
    

    #***************************************************************************************
    # Estimate initial pmf
    #***************************************************************************************
    #the initial pmf should be DGATP/N
    n=4.666
    pmf_initial=0.0
    #***************************************************************************************
    #Initial lumen pH
    #the following sets up the initial pmf and pH values
    #pH_stroma will be held constant
    #***************************************************************************************

    pHstroma_initial=7.8
    #pHstroma=pHstroma_initial

    pHlumen_initial=7.8 #initially, place abouit half of pmf into DpH

    #***************************************************************************************
    #Initial Dy
    #the following sets up the initial Dy
    #***************************************************************************************

    Dy_initial=0.0 #place the other half as Dy

    LEF_initial=0
    Phi2_initial=0.83


    #print('With n=' + str(n) + ', the pmf at equilibrium with DGATP should be set to : ' + str(pmf_initial))

    #tell the user what the concentration of free H+ is in the lumen
    #free_H=10**(-1*pmf_initial)
    #print('the estimated concentration of free protons in the lumen = '  + str(free_H))

    #tell the user the concentration of total protons in the lumen
    buffering_capacity=0.03
    Hin_initial= 0.0
    Hstroma_initial = 0.0
    #print('the concentration of total (free + bound) protons in the lumen = ' + str(Hin_initial))

    #pHlumen_initial=7-Hin_initial/buffering_capacity
    #print('the initial lumen pH = ' + str(pHlumen_initial))
    QA_content_initial=1
    QAm_content_initial=0

    #***************************************************************************************
    #parameters for PQ pool 
    #***************************************************************************************

    PQH2_content_initial=0
    PQ_content_initial=7

    #***************************************************************************************
    #parameters for Plastocyanin (PC) reactions 
    #***************************************************************************************

    PC_ox_initial = 0
    PC_red_initial= 2

    #***************************************************************************************
    #parameters for P700 reactions 
    #***************************************************************************************

    P700_ox_initial=0.0
    P700_red_initial=0.667 #Mathias Pribil1, Mathias Labs1 and Dario Leister1,2,* Structure and dynamics of thylakoids in land plantsJournal of Experimental Botany, Vol. 65, No. 8, pp. 1955–1972, 2014
    #Fan DY1, Hope AB, Smith PJ, Jia H, Pace RJ, Anderson JM, Chow WS, The stoichiometry of the two photosystems in higher plants revisited. Biochim Biophys Acta. 2007 Aug;1767(8):1064-72

    PSI_content=P700_red_initial + P700_ox_initial#PSI/PSII = 0.75, see doi:10.1093/jxb/eru090
    PSI_antenna_size=0.5 #setting this to the same valus as PSII_antenna_size will imply that 
                        #equal numbers of photons hit the two photosystems


    Fd_ox_initial=1 #currently, Fd will just accumulate 
    Fd_red_initial=0 #currently, Fd will just accumulate 
    

    #***************************************************************************************
    #parameters for NPQ 
    #***************************************************************************************
    NPQ_initial=0

    singletO2_initial=0
    ATP_pool_initial=4.15#dark ATP/ADP = 1, light ~5, 1 mM ATP under light 1.5 mM Pi constant
    ADP_pool_initial=4.15
    NADPH_pool_initial=1.5
    NADP_pool_initial=3.5
    
#*******************************************************************************
#*******************************************************************************
#                    Classes to hold constants and states. 
#*******************************************************************************
#********************************************************************************

class sim_constants(FrozenClass):
    def __init__(self):
        S=standard_constants()
        self.pKreg=S.pKreg
        self.max_PSII=S.max_PSII
        self.kQA=S.kQA
        self.max_b6f=S.max_b6f
        self.lumen_protons_per_turnover=S.lumen_protons_per_turnover
        self.light_per_L=S.light_per_L
        self.ATP_synthase_max_turnover=S.ATP_synthase_max_turnover
        #self.pHstroma=S.pHstroma_initial
        self.PSII_antenna_size=S.PSII_antenna_size
        self.Volts_per_charge=S.Volts_per_charge
        self.perm_K=S.perm_K
        self.n=S.n
        self.Em7_PQH2=S.Em7_PQH2
        self.Em7_PC=S.Em7_PC
        self.Em_Fd = S.Em_Fd
        self.PSI_antenna_size=S.PSI_antenna_size
        self.buffering_capacity=S.buffering_capacity
        self.VDE_max_turnover_number=S.VDE_max_turnover_number
        self.pKvde=S.pKvde
        self.VDE_Hill=S.VDE_Hill
        self.kZE=S.kZE
        self.pKPsbS=S.pKPsbS
        self.max_NPQ=S.max_NPQ
        self.k_recomb=S.k_recomb
        self.k_PC_to_P700=S.k_PC_to_P700
        self.triplet_yield=S.triplet_yield
        self.triplet_to_singletO2_yield=S.triplet_to_singletO2_yield
        self.fraction_pH_effect=0.25
        self.k_Fd_to_NADP=S.k_Fd_to_NADP
        self.k_CBC=S.k_CBC
        self.k_KEA=S.k_KEA
        self.k_VCCN1 = S.k_VCCN1
        self.k_CLCE = S.k_CLCE
        self.k_NDH = S.k_NDH
        
    def as_tuple(self):
        c=(self.pKreg, self.max_PSII, self.kQA, self.max_b6f, self.lumen_protons_per_turnover, self.light_per_L,
        self.ATP_synthase_max_turnover, self.PSII_antenna_size, self.Volts_per_charge, self.perm_K,
        self.n, self.Em7_PQH2, self.Em7_PC, self.Em_Fd, self.PSI_antenna_size, self.buffering_capacity, 
        self.VDE_max_turnover_number, self.pKvde, self.VDE_Hill, self.kZE, self.pKPsbS, self.max_NPQ, 
        self.k_recomb, self.k_PC_to_P700, self.triplet_yield, self.triplet_to_singletO2_yield, 
        self.fraction_pH_effect, self.k_Fd_to_NADP, self.k_CBC, self.k_KEA, self.k_VCCN1, self.k_CLCE, self.k_NDH)
        return(c)
    
    def as_dictionary(self):
        
        d={'pKreg':self.pKreg, 'max_PSII':self.max_PSII,'kQA': self.kQA, 'max_b6f': self.max_b6f, 
        'lumen_protons_per_turnover': self.lumen_protons_per_turnover, 'light_per_L':self.light_per_L,
        'ATP_synthase_max_turnover': self.ATP_synthase_max_turnover,  
        'PSII_antenna_size': self.PSII_antenna_size, 'Volts_per_chargese': self.Volts_per_charge, 
        'perm_K': self.perm_K, 'n': self.n, 'Em7_PQH2': self.Em7_PQH2, 'Em7_PC': self.Em7_PC,'Em_Fd':self.Em_Fd, 
        'PSI_antenna_size': self.PSI_antenna_size, 'buffering_capacity': self.buffering_capacity, 
        'VDE_max_turnover_number': self.VDE_max_turnover_number, 'pKvde': self.pKvde, 'VDE_Hill': self.VDE_Hill, 
        'kZE': self.kZE, 'pKPsbS': self.pKPsbS, 'max_NPQ': self.max_NPQ, 
        'k_recomb': self.k_recomb, 'k_PC_to_P700': self.k_PC_to_P700, 'triplet_yield': self.triplet_yield, 
        'triplet_to_singletO2_yield': self.triplet_to_singletO2_yield, 'fraction_pH_effect': self.fraction_pH_effect, 
        "k_Fd_to_NADP":self.k_Fd_to_NADP, "k_CBC": self.k_CBC, "k_KEA":self.k_KEA, 'k_VCCN1':self.k_VCCN1,
        'k_CLCE':self.k_CLCE, 'k_NDH':self.k_NDH}
        return(d)
        
    def short_descriptions(self):
        e={'pKreg':'The regulatory pKa at which the cytochrome b6f complex is slowed by lumen pH', 
        'max_PSII':'The maximum relative rate of PSII centers (0-1)',
        'kQA': 'The averate rate constant of oxidation of QA- by QB and PQ', 
        'max_b6f': 'The maximum turnover rate for oxidation of PQH2 by the b6f complex at high pH', 
        'lumen_protons_per_turnover': 'The molarity change of protons resulting from 1 H+/standard PSII into the lumen', 
        'light_per_L':'PAR photons per PSII',
        'ATP_synthase_max_turnover': 'Defines the slope of ATP synthesis per pmf', 
        
        'PSII_antenna_size': 'The relative antenna size of PSII', 
        'Volts_per_chargese': 'The capcitance of the thylakoid expressed as V/charge/PSII', 
        'perm_K': 'The relative permeability of the thylakoid to counterions', 
        'n': 'The stoichiometry of H+/ATP at the ATP synthase', 
        'Em7_PQH2': 'The midpoint potential of the PQ/PQH2 couple at pH=7', 
        'Em7_PC': 'The midpoint potential of the plastocyanin couple at pH=7',
        'Em_Fd':'The midpoint potential of ferredoxin',
        'PSI_antenna_size': 'The relative cross section of PSI antenna', 
        'buffering_capacity': 'The proton buffering capacity of the lumen in M/pH unit', 
        'VDE_max_turnover_number': 'The maximal turnover of the fully protonated VDE enzyme', 
        'pKvde': 'The pKa for protonation and activation of VDE', 
        'VDE_Hill': 'The Hill coefficient for protonation of VDE', 
        'kZE': 'The rate constant for ZE (zeaxanthin epoxidase', 
        'pKPsbS': 'The pKa for protonation and activation of PsbS', 
        'max_NPQ': 'NPQ=(nax_NPQ)(PsbS protonation)(Z)', 
        'k_recomb': 'The average recombination rate for S2QA- and S3QA- with no detla.psi (s-1)', 
        'k_PC_to_P700': 'The rate constant for transfer of electrons from PC to P700+',
        'triplet_yield': 'The yield of triplets from each recombination event', 
        'triplet_to_singletO2_yield': 'The yield of 1O2 for each triplet formed', 
        'fraction_pH_effect': 'The frqaction of S-states that both involve protons release and can reconbine', 
        "k_Fd_to_NADP":'The second order rate constant for oxidation of Fd by NADP+', 
        "k_CBC": 'The rate constant for a simplified Calvin-Benson Cycle',
        "k_KEA": 'The rate constant for the KEA H+/H+ antiporter',
        'k_VCCN1':'The rate constant for VCCN1 moving Cl- from stroma to lumen',
        'k_CLCE':'The rate constant for CLCE2 moving Cl- from lumen to stroma, driving by H+ gradient',
        'k_NDH':'The rate constant for NDH'}
        return(e)
        
        #self._freeze() # no new attributes after this point.
        
        
class sim_states(FrozenClass):
    def __init__(self):
        S=standard_initial_states()
        self.QA_content=S.QA_content_initial
        self.QAm_content=S.QAm_content_initial
        self.PQ_content=S.PQ_content_initial
        self.PQH2_content=S.PQH2_content_initial
        self.Hin=S.Hin_initial
        self.pHlumen=S.pHlumen_initial
        self.Dy=S.Dy_initial
        self.pmf=S.pmf_initial
        self.DeltaGatp=S.DeltaGatp_initial
        self.Klumen=S.Klumen_initial
        self.Kstroma=S.Kstroma_initial
        self.ATP_made=S.ATP_made_initial
        self.PC_ox=S.PC_ox_initial
        self.PC_red=S.PC_red_initial
        self.P700_ox=S.P700_ox_initial
        self.P700_red=S.P700_red_initial
        self.Z=S.Z_initial
        self.V=S.V_initial
        self.NPQ=S.NPQ_initial
        self.singletO2=S.singletO2_initial
        self.Phi2=S.Phi2_initial
        self.LEF=S.LEF_initial
        self.Fd_ox=S.Fd_ox_initial
        self.Fd_red=S.Fd_red_initial
        self.ATP_pool=S.ATP_pool_initial
        self.ADP_pool=S.ADP_pool_initial
        self.NADPH_pool=S.NADPH_pool_initial
        self.NADP_pool=S.NADP_pool_initial
        self.Cl_lumen = S.Cl_lumen_initial
        self.Cl_stroma = S.Cl_stroma_initial
        self.H_stroma = S.Hstroma_initial
        self.pHstroma = S.pHstroma_initial
        
    def as_list(self):
            t=[self.QA_content, self.QAm_content, self.PQ_content, 
                       self.PQH2_content, self.Hin, self.pHlumen, self.Dy, self.pmf, self.DeltaGatp,
                       self.Klumen, self.Kstroma, self.ATP_made, self.PC_ox, 
                       self.PC_red, self.P700_ox, self.P700_red, self.Z,self.V, self.NPQ,
                       self.singletO2, self.Phi2, self.LEF, self.Fd_ox, self.Fd_red, self.ATP_pool, 
                       self.ADP_pool, self.NADPH_pool, self.NADP_pool, self.Cl_lumen, self.Cl_stroma,
                       self.H_stroma, self.pHstroma]
            return(t)
        
    def as_tuple(self):
            t=tuple([self.QA_content, self.QAm_content, self.PQ_content, 
                       self.PQH2_content, self.Hin, self.pHlumen, self.Dy, self.pmf, self.DeltaGatp,
                       self.Klumen, self.Kstroma, self.ATP_made, self.PC_ox, 
                       self.PC_red, self.P700_ox, self.P700_red, self.Z,self.V, self.NPQ,
                       self.singletO2, self.Phi2, self.LEF, self.Fd_ox, self.Fd_red, self.ATP_pool, 
                       self.ADP_pool, self.NADPH_pool, self.NADP_pool, self.Cl_lumen, self.Cl_stroma,
                       self.H_stroma, self.pHstroma])
            return(t)
        

    def load_from_tupple(self, Y):
        self.QA_content=Y[0]
        self.QAm_content=Y[1]
        self.PQ_content=Y[2] 
        self.PQH2_content=Y[3]
        self.Hin=Y[4]
        self.pHlumen=Y[5]
        self.Dy=Y[6]
        self.pmf=Y[7]
        self.DeltaGatp=Y[8]
        self.Klumen=Y[9]
        self.Kstroma=Y[10]
        self.ATP_made=Y[11]
        self.PC_ox=Y[12]
        self.PC_red=Y[13]
        self.P700_ox=Y[14]
        self.P700_red=Y[15]
        self.Z=Y[16]
        self.V=Y[17]
        self.NPQ=Y[18]
        self.singletO2=Y[19]
        self.Phi2=Y[20]
        self.LEF=Y[21]
        self.Fd_ox=Y[22]
        self.Fd_red=Y[23]
        self.ATP_pool=Y[24]
        self.ADP_pool=Y[25]
        self.NADPH_pool=Y[26]
        self.NADP_pool=Y[27]
        self.Cl_lumen = Y[28]
        self.Cl_stroma = Y[29]
        self.H_stroma = Y[30]
        self.pHstroma = Y[31]
    
    def as_dictionary(self):
        d={'QA_content': self.QA_content,
        'QAm_content':self.QAm_content,
        'PQ_content':self.PQ_content,
        'PQH2_content': self.PQH2_content,
        'Hin':self.Hin,
        'pHlumen': self.pHlumen,
        'Dy':self.Dy,
        'pmf':self.pmf,
        'DeltaGatp': self.DeltaGatp,
        'Klumen':self.Klumen,
        'Kstroma': self.Kstroma,
        'ATP_made': self.ATP_made,
        'PC_ox': self.PC_ox,
        'PC_red': self.PC_red,
        'P700_ox':self.P700_ox,
        'P700_red': self.P700_red,
        'Z':self.Z,
        'V':self.V,
        'NPQ':self.NPQ,
        'singletO2':self.singletO2,
        'Phi2':self.Phi2,
        'LEF': self.LEF,
        'Fd_ox': self.Fd_ox,
        'Fd_red': self.Fd_red,
        'ATP_pool': self.ATP_pool,
        'ADP_pool':self.ADP_pool,
        'NADPH_pool': self.NADPH_pool,
        'NADP_pool':self.NADP_pool,
        'Cl_lumen':self.Cl_lumen,
        'Cl_stroma':self.Cl_stroma,
        'Hstroma':self.H_stroma,
        'pHstroma':self.pHstroma}
        self._freeze() # no new attributes after this point.
        return(d)

    
#*******************************************************************************
#*******************************************************************************
#                    Display notes, states and constants. 
#*******************************************************************************
#********************************************************************************


class ListTable(list):
    """ Overridden list class which takes a 2-dimensional list of 
        the form [[1,2,3],[4,5,6]], and renders an HTML Table in 
        IPython Notebook. """
    
    def _repr_html_(self):
        html = ["<table>"]
        for row in self:
            html.append("<tr>")
            
            for col in row:
                html.append("<td>{0}</td>".format(col))
            
            html.append("</tr>")
        html.append("</table>")
        return ''.join(html)
        
#Display all constants in a table
def All_Constants_Table(table_title, Kxx):
    table = ListTable()
    table.append(['Parameter', 'New Value']) #, 'Short Description'])
    #Kxx=sim_constants()
    Kdict=Kxx.as_dictionary()
    #Ksumm=Kxx.short_descriptions()

    for key in list(Kdict.keys()):

        table.append([key, Kdict[key]]) #, Ksumm[key]])
    print(table_title)
    display(table)
    
#display only the constants that are different
def Changed_Constants_Table(table_title, original_values, Kxx):
    table = ListTable()
    table.append(['Changed Parameter', 'Old Value', 'New Value']) #, 'Short Description'])
    Kdict=Kxx.as_dictionary()
    #Ksumm=Kxx.short_descriptions()
    #temp=sim_constants()
    original_values_dict=original_values.as_dictionary()

    for key in list(Kdict.keys()):
        if Kdict[key] == original_values_dict[key]:
            pass
        else:
            table.append([key, original_values_dict[key], Kdict[key]]) #, Ksumm[key]])
    print(table_title)
    display(table)



#Display PDFs of the detailed notes describing the simulation parameters
#def display_detailed_notes():
#    from IPython.display import IFrame
#    from IPython.display import Image
#
#    page1=Image(PDF_file_location + 'Page 1.png')
#    page2=Image(PDF_file_location + 'Page 2.png')
#    page3=Image(PDF_file_location + 'Page 3.png')

#    display(page1)
#    display(page2)
#    display(page3)

#SAVE A LIST(SIMULATED VALUES) TO CSV FILE
def list_to_csv(alist, filename):
    with open(filename,'w',newline='') as f:
        writer_a = csv.writer(f)
        for aline in alist:
            writer_a.writerow(aline)
#PROCESS A GENOTYPE AND SAVE ITS SIMULATED VALUES
def process_a_gtype(gtype_dict, parameter_list, out_dict, gtype='a_genotype'):
    gtype_df = pd.DataFrame([])
    for para in parameter_list:
        gtype_dict[para] = out_dict[para]#store in dictionary for further calculation
        # if type(gtype_dict[para]) is list:
        #     temp_list0 = gtype_dict[para]
        # else:
        #     temp_list0 = gtype_dict[para].tolist()
        # temp_list = [para]+temp_list0[:]
        #this slice is due to too many data points at transition from dark to light, & light to dark
        gtype_df[para] = out_dict[para]
    #list_to_csv(gtype_list, gtype+'_simulated.csv')
    gtype_df.to_csv(gtype+'_simulated.csv')

#the following function simulate a specific genotype, change ATPase activation function   
def sim_a_gtype(gtype_dict, gtype='WT', light = 100):  
    parameters_of_interest = ['time_axis','NPQ','Phi2','LEF','qL','Z','V',\
                          'pmf','Dy','pHlumen','fraction_Dy','fraction_DpH',\
                          'Klumen','Cl_lumen','Cl_stroma']
    #this parameters_of_interest is a list of things exported into csv and compared
    #between wt and mutants, do_complete_sim dictates how each paarmeter is called
    #run the code to make all pre-contrived light waves
    #light_pattern=make_waves()    
    #The following code generates a series of diurnal light patters, with either smooth or fluctuating
    #patterns
    initial_sim_states=sim_states()
    initial_sim_state_list=initial_sim_states.as_list()
    Kx_initial=sim_constants()    
    #All_Constants_Table('Standard Constants', Kx_initial)
    constants_dict={}
    #starting_conditions_dict={}
    k_CBC_light = 60 * (light/(light+250))#this needs change with different light intensity    
    ####this following name WT as a dictionary, when WT[parameter] is called,
    ####it will return the parameter np_array, and is convenient for calculations
    output_dict={}
    on = gtype
    Kx=sim_constants()
    if 'clce2' in gtype:
        Kx.k_CLCE = 0
    if 'kea3' in gtype:
        Kx.k_KEA =0
    if 'vccn1' in gtype:
        Kx.k_VCCN1 =0
    Kx.k_CBC = k_CBC_light
    constants_dict[on]=Kx #store constants in constants_dict
    # if light ==100:
    #     output_dict, starting_conditions_dict[on]=sim(Kx, initial_sim_state_list, 
    #                                         light_pattern['single_square_20_min_100_max'])
    #if light == 500:
    output_dict=sim_ivp(Kx, initial_sim_state_list, 1200)
    Changed_Constants_Table('Change Constants', Kx_initial, Kx)
    output_dict['qL'] = 1-output_dict['QAm']
    plot_interesting_stuff(gtype, output_dict)
    process_a_gtype(gtype_dict,parameters_of_interest, output_dict,gtype+'_'+str(light)+'uE')    
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


    


