import numpy as np
import pandas as pd


class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

        
    def as_tuple(self):
        c=(self.ratio_absorb,self.PSII_content_per_leaf,self.pKreg, self.max_PSII, self.kQA, self.max_b6f, self.lumen_protons_per_turnover, self.light_per_L,
        self.ATP_synthase_max_turnover, self.PSII_antenna_size, self.Volts_per_charge, self.perm_K,
        self.n, self.Em7_PQH2, self.Em7_PC, self.Em_Fd, self.PSI_antenna_size, self.buffering_capacity, 
        self.VDE_max_turnover_number, self.pKvde, self.VDE_Hill, self.kZE, self.pKPsbS, self.max_NPQ, 
        self.k_recomb, self.k_PC_to_P700, self.triplet_yield, self.triplet_to_singletO2_yield, 
        self.fraction_pH_effect, self.k_Fd_to_NADP, self.k_CBC, self.k_KEA, self.k_VCCN1, self.k_CLCE, self.k_NDH)
        return(c)
    
    def as_dictionary(self):
        d={'ratio_absorb':self.ratio_absorb,'PSII_content_per_leaf':self.PSII_content_per_leaf,'pKreg':self.pKreg, 'max_PSII':self.max_PSII,'kQA': self.kQA, 'max_b6f': self.max_b6f, 
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

class standard_constants(object):
    def __init__(self, csv_file='./data/constants.csv'):
        if csv_file:
            self.load_from_csv(csv_file)
        else:
            print("Reading CSV file error, check CSV file!!!")
    def load_from_csv(self, csv_file):
        df = pd.read_csv(csv_file)
        for column in df.columns:
            setattr(self, column, df[column].iloc[0])

class standard_initial_states(object):
    def __init__(self, csv_file='./data/initial_states.csv'):
        if csv_file:
            self.load_from_csv(csv_file)
        else:
            print("Reading CSV file error, check CSV file!!!")

    def load_from_csv(self, csv_file):
        df = pd.read_csv(csv_file)
        for column in df.columns:
            setattr(self, column, df[column].iloc[0])

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

class sim_constants(FrozenClass):
    def __init__(self):
        S=standard_constants()
        self.ratio_absorb=S.ratio_absorb
        self.PSII_content_per_leaf=S.PSII_content_per_leaf
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