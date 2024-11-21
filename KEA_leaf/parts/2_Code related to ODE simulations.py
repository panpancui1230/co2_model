# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 19:43:01 2024

@author: LENOVO
"""

   

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

