import numpy as np

import componentmanager
from flow import *
from componentmanager import GlobalComponentManager
from reactions import *
from proptiesmethod import *


class CstrSingleLiqPhase:

    def __init__(self, t, p, v, inflow: Flow, rx_set: ReactionSet, prop: PropertiesMethod):

        self.ReactionSet = rx_set
        self.Inflow = inflow
        self.PropertiesMethod = prop
        self.Temperature = t
        self.Pressure = p
        self.Volume = v
        self.Mole_Frac_Dict = {}
        self.mole_frac_init(3)

    def mole_frac_init(self, site):

        self.Mole_Frac_Dict = {c: 0. for c in GlobalComponentManager.component_list}
        self.Mole_Frac_Dict['ion-pair'] = 0.
        active_spices = []
        p1_ion = []

        zeroth_mom_live = []
        first_mom_live = []
        second_mom_live = []

        zeroth_mom_dead = []
        first_mom_dead = []
        second_mom_dead = []

        for i in range(0, site):
            active_spices.append(0.0)
            p1_ion.append(0.0)
            zeroth_mom_live.append(0.0)
            first_mom_live.append(0.0)
            second_mom_live.append(0.0)

            zeroth_mom_dead.append(0.0)
            first_mom_dead.append(0.0)
            second_mom_dead.append(0.0)

        self.Mole_Frac_Dict['active_species'] = active_spices
        self.Mole_Frac_Dict['p1_ion'] = p1_ion
        self.Mole_Frac_Dict['zeroth_mom_live'] = zeroth_mom_live
        self.Mole_Frac_Dict['first_mom_live'] = first_mom_live
        self.Mole_Frac_Dict['second_mom_live'] = second_mom_live
        self.Mole_Frac_Dict['zeroth_mom_dead'] = zeroth_mom_dead
        self.Mole_Frac_Dict['first_mom_dead'] = first_mom_dead
        self.Mole_Frac_Dict['second_mom_dead'] = second_mom_dead
        self.Mole_Frac_Dict['counter_ion'] = 0.0

    def mass_balance(self, Outflow: Flow):

        mb_eq = []

        polymer_Mole_flow_zeroth = np.array(
            [mole['mole_flow'] if mole['polymer_flow_momentum'] == 0 else 0.0
             for mole in Outflow.comp_dict])

        polymer_Mole_flow_first = np.array(
            [mole['mole_flow'] if mole['polymer_flow_momentum'] == 1 else 0.0
             for mole in Outflow.comp_dict])

        Mole_flow_zeroth = np.array(
            [mole['mole_flow'] if mole['polymer_flow_momentum'] in [0, -2] else 0.0 for
             mole in Outflow.comp_dict])

        Mole_flow_first = np.array(
            [mole['mole_flow'] if mole['polymer_flow_momentum'] in [1, -2] else 0.0 for
             mole in Outflow.comp_dict])

        Mole_frac_zeroth = Mole_flow_zeroth / np.sum(Mole_flow_zeroth)
        Mole_frac_first = Mole_flow_first / np.sum(Mole_flow_first)
        dpn = (np.sum(polymer_Mole_flow_first) / np.sum(polymer_Mole_flow_zeroth)) if np.sum(
            polymer_Mole_flow_zeroth) != 0. else 1e12

        for c_idx, comp in enumerate(GlobalComponentManager.component_list):
            if type(comp) is Component:
                pass

        self.PropertiesMethod.param[-1]=0.

        self.PropertiesMethod.calculate_molar_density_mixture(self.Temperature, self.Pressure)
        # molar density [mol/m3]
        vm_liq = self.PropertiesMethod.calculate_molar_density_mixture(self.Mole_Frac_Dict, self.Temperature,
                                                                       self.Pressure)

        q_out = np.sum(Mole_flow_first) / vm_liq

        concentration = np.array([c['mole_flow']/q_out for c in Outflow.comp_dict])

        for f in self.Inflow.comp_dict:
            eq = self.Inflow.comp_dict[f]['mole_flow'] - Outflow.comp_dict[f]['mole_flow'] + self.ReactionSet.calculate_rate(f,
                                                                                                   concentration) * self.Volume
            mb_eq.append(eq)

        return mb_eq
