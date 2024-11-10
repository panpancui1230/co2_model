import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter

class sim_plot:
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
        self._freeze()

class Plotting:
    def __init__(self):
        pass

    def plot_interesting_stuff(self, figure_name, output):
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

    def plot_gen(self, fig, sub_plot_number, plot_list, plot_every_nth_point, **keyword_parameters):
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

    def plot_pmf_params(self, output, use_x_axis, x_axis_label, global_min, global_max):
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

    def plot_pmf_params_offset(self, output, use_x_axis, x_axis_label, global_min, global_max):
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

    def plot_K_and_parsing(self, output, use_x_axis, x_axis_label, global_min, global_max):
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

    def b6f_and_balance(self, output, use_x_axis, x_axis_label, global_min, global_max):
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

    def plot_QAm_and_singletO2(self, output, use_x_axis, x_axis_label, global_min, global_max):
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

    def plot_cum_LEF_singletO2(self, output, use_x_axis, x_axis_label, global_min, global_max):
        # set up the left axis of the plot for commulative LEF
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

    def plot_differences(self, time_min, delta_NPQ, delta_LEF):
        fig = plt.figure(num=3, figsize=(5, 4), dpi=200)
        plt.plot(time_min[1:], delta_NPQ[1:], label='∆NPQ: kea3 - WT')
        plt.legend()
        plt.show()
        plt.close()

        fig = plt.figure(num=3, figsize=(5, 4), dpi=200)
        plt.plot(time_min[1:], delta_LEF[1:], label='∆LEF: kea3 - WT')
        plt.legend()
        plt.show()
        plt.close()

    def plot_cum_LEF_singetO2(self, output, use_x_axis, x_axis_label,global_min, global_max):
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

    def b6f_and_balance(self, output, use_x_axis, x_axis_label, global_min, global_max):

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