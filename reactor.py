import componentmanager
from flow import *
from componentmanager import GlobalComponentManager
from reactions import *
from poptiesmethod import *

class cstr_SingleLiqPhase():

    def __init__(self, T, p, V, Inflow:Flow,rx_set:ReactionSet,prop:PropertiesMethod):
        self.ReactionSet=rx_set
        self.Inflow=Inflow
        self.PropertiesMethod=prop
        self.Temperature = T
        self.Pressure = p
        self.Volume = V
        self.Mole_Frac_Dict = {}
        self.mole_frac_init(3)

    def mole_frac_init(self,site):

        self.Mole_Frac_Dict = {c: 0. for c in GlobalComponentManager.component_list}
        self.Mole_Frac_Dict['ion-pair']=0.
        active_spices=[]
        p1_ion=[]

        zeroth_mom_live=[]
        first_mom_live=[]
        second_mom_live=[]

        zeroth_mom_dead=[]
        first_mom_dead=[]
        second_mom_dead=[]

        for i in range(0,site):
            active_spices.append(0.0)
            p1_ion.append(0.0)
            zeroth_mom_live.append(0.0)
            first_mom_live.append(0.0)
            second_mom_live.append(0.0)

            zeroth_mom_dead.append(0.0)
            first_mom_dead.append(0.0)
            second_mom_dead.append(0.0)

        self.Mole_Frac_Dict['active_species']=active_spices
        self.Mole_Frac_Dict['p1_ion'] = p1_ion
        self.Mole_Frac_Dict['zeroth_mom_live']=zeroth_mom_live
        self.Mole_Frac_Dict['first_mom_live']=first_mom_live
        self.Mole_Frac_Dict['second_mom_live']=second_mom_live
        self.Mole_Frac_Dict['zeroth_mom_dead']=zeroth_mom_dead
        self.Mole_Frac_Dict['first_mom_dead']=first_mom_dead
        self.Mole_Frac_Dict['second_mom_dead']=second_mom_dead
        self.Mole_Frac_Dict['counter_ion']=0.0


    def mass_balance(self,solu):

        MB_eq=[]
        vm=self.PropertiesMethod.calculate_molar_density_mixture(self.Temperature)

        if vm!=0:
            concen = [s / vm for s in solu]
        else:
            concen=[0. for i in range(len(solu))]

        for f in self.Inflow.comp_dict:
            eq= self.Inflow.comp_dict[f]['mole_flow']-solu[0]+self.ReactionSet.calculate_rate(f,concen)*self.Volume
            MB_eq.append(eq)

        return MB_eq

