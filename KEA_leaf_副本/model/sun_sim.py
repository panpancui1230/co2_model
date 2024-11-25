import numpy as np

class sunshine:
    def optimized_time_split(self, test_times_and_light, max_light_change, points_per_segment):   
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

    def make_variable_light_constants_set_and_trace_times(self, K, sub_arrays_time_and_light):
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

    def generate_sin_wave(self, total_duration, max_PAR, light_frequency, points_per_second):
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
    
    def square_light(self, t, light_intensity, duration = 20, unit= 'min',\
                 t0 = 0, prior_light = 0, post_light = 0 ):

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
    
    def sin_light_fluctuating(self, t, freq, PAR_max, PAR_min, t0=0, PAR_0=None ):
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

    def light(self, t, Duration_T0, par0, frequency, par_max, par_min):
        if t <= Duration_T0:
            par = self.square_light(t, par0, Duration_T0, 'seconds')
        else:
            par = self.sin_light_fluctuating(t, frequency, par_max, par_min, \
                                        t0 = Duration_T0, PAR_0 = par0)
        return par