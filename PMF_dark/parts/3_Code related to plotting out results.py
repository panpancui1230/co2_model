# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:45:43 2024

@author: LENOVO
"""


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
    
