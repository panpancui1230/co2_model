# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:47:57 2024

@author: LENOVO
"""


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
