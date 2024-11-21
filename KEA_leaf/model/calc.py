import numpy as np

class block:
    def calc_k_b6f(self, max_b6f, b6f_content, pHlumen, pKreg):
        pHmod=(1 - (1 / (10 ** (pHlumen - pKreg) + 1)))
        b6f_deprot=pHmod*b6f_content
        k_b6f=b6f_deprot * max_b6f
        return(k_b6f)
    
    def calc_v_b6f(self, max_b6f, b6f_content, pHlumen, pKreg, PQ, PQH2, PC_ox, PC_red, Em7_PC, Em7_PQH2, pmf):    
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
    
    def calc_v_NDH(self, Em_Fd, Em7_PQH2, pHstroma, pmf, k_NDH, Fd_red, Fd_ox, PQ, PQH2):
        Em_PQH2 = Em7_PQH2 - 0.06*(pHstroma - 7.0)
        deltaEm = Em_PQH2 - Em_Fd
        Keq_NDH = 10**((deltaEm - pmf*2)/0.06)
        k_NDH_reverse = k_NDH/Keq_NDH
        #f_PQ = PQ/(PQ+PQH2)
        #f_PQH2 = 1.0-f_PQ
        v_NDH = k_NDH*Fd_red*PQ - k_NDH_reverse*Fd_ox*PQH2
        return (v_NDH)
    
    def calc_v_PGR(self, PGR_vmax, Fd_red, PQ, PQH2):
        v_PGR = PGR_vmax * (Fd_red**4/(Fd_red**4+0.1**4))*PQ/(PQ+PQH2)
        return v_PGR
    
    #calculate the rate of V<-- -->Z reactions, assuming a pH-dependent VDE and a pH-independent ZE
    def calc_v_VDE(self, VDE_max_turnover_number, pKvde, VDE_Hill, kZE, pHlumen, V, Z):    
        pHmod= 1 / (10 ** (VDE_Hill*(pHlumen - pKvde)) + 1)
        v_Z = V* VDE_max_turnover_number*pHmod - Z*kZE
        v_V = -1* v_Z
        return(v_Z, v_V)
    
    def calc_CBC_NADPH(self, k_CBC,t,d_ATP_made):
        NADPH_CBC_t =  k_CBC*(1.0-np.exp(-t/600))
        NADPH_CBC_ATP =0.6*d_ATP_made
        NADPH_CBC = min([NADPH_CBC_ATP,NADPH_CBC_t])
        return NADPH_CBC
            
    def Calc_Phi2(self, QA, NPQ):
        Phi2=1/(1+(1+NPQ)/(4.88*QA))
        return Phi2

    def Calc_PhiNO_PhiNPQ(self, Phi2, QA, NPQ):
        PhiNO=1/(1+NPQ + ((Phi2+NPQ)/(1-Phi2)))
        PhiNPQ=1-(Phi2+PhiNO)
        return PhiNO, PhiNPQ
    
    def calc_PsbS_Protonation(self, pKPsbS, pHlumen):    
        PsbS_H=1 / (10 ** (3*(pHlumen - pKPsbS)) + 1)
        return(PsbS_H)
    
    def ATP_synthase_actvt(self, t, T_ATP):#based on gH+ data
        x = t/T_ATP
        actvt = 0.2 + 0.8*(x**4/(x**4 + 1))
        return actvt
    
    def Vproton_pmf_actvt(self, pmf, actvt, ATP_synthase_max_turnover, n):# fraction of activity based on pmf, pmf_act is the half_max actvt pmf
        v_proton_active = 1 - (1 / (10 ** ((pmf - 0.132)*1.5/0.06) + 1))#reduced ATP synthase
        v_proton_inert = 1-(1 / (10 ** ((pmf - 0.204)*1.5/0.06) + 1))#oxidized ATP synthase
        
        v_active = actvt * v_proton_active * n * ATP_synthase_max_turnover
        v_inert = (1-actvt) * v_proton_inert * n * ATP_synthase_max_turnover
        
        v_proton_ATP = v_active + v_inert
        return (v_proton_ATP)

    def V_H_light(self, light_per_L, v_proton_ATPase, pmf, Hlumen, k_leak = 3*10**7):
        if light_per_L>0.0:
            V_H = v_proton_ATPase + pmf*k_leak*Hlumen
        else:
            V_H = pmf*k_leak*Hlumen# this term is used for dark relaxation,
            #ATP synthase actvt dependent but does not make ATP
        return V_H

    def recombination_with_pH_effects(self, k_recomb, QAm, Dy, pHlumen, fraction_pH_effect):
        delta_delta_g_recomb= Dy + .06*(7.0-pHlumen)
        v_recomb = k_recomb*QAm*10**(delta_delta_g_recomb/.06)
        
        #v_recomb = k_recomb*QAm*(10**((Dy/.06) + fraction_pH_effect*10**(7.0-pHlumen)))        
        return(v_recomb)

    def Cl_flux_relative(self, v):
        Cl_flux_v = 332*(v**3) + 30.8*(v**2) + 3.6*v
        #relative to Cl flux thru VCCN1. when driving force is 0.1 Volt,
        #Cl_flux_v is 1. empirical equation was obtained from
        # Herdean et al. 2016 DOI: 10.1038/ncomms11654
        return Cl_flux_v
    def KEA_reg(self, pHlumen, QAm):
        qL = 1-QAm
        qL_act = qL**3/(qL**3+0.15**3)
        pH_act =1/(10**(1*(pHlumen-6.0))+1)
        f_KEA_act = qL_act * pH_act
        return f_KEA_act