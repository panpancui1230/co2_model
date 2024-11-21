# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:52:50 2024

@author: LENOVO
"""

    
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
    
 

    


