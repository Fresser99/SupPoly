import numpy as np

import componentmanager
from flow import *
from componentmanager import GlobalComponentManager
from reactions import *
from proptiesmethod import *
import pyomo.environ as pyo


class CstrSingleLiqPhase:

    def __init__(self, t, p, v, inflow: Flow, rx_set: ReactionSet, prop: PropertiesMethod):

        self.ReactionSet = rx_set
        self.Inflow = inflow
        self.PropertiesMethod = prop
        self.Temperature = t
        self.Pressure = p
        self.Volume = v

    def mass_balance(self, Outflow: Flow):

        mb_eq = []

        polymer_Mole_flow_zeroth = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole]['polymer_flow_momentum'] == 0 else 0.0
             for mole in Outflow.comp_dict])

        polymer_Mole_flow_first = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole]['polymer_flow_momentum'] == 1 else 0.0
             for mole in Outflow.comp_dict])

        Mole_flow_zeroth = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole]['polymer_flow_momentum'] in [0,
                                                                                                          -2] else 0.0
             for
             mole in Outflow.comp_dict])

        Mole_flow_first = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole]['polymer_flow_momentum'] in [1,
                                                                                                          -2] else 0.0
             for
             mole in Outflow.comp_dict])

        Mole_frac_zeroth = Mole_flow_zeroth / np.sum(Mole_flow_zeroth)
        Mole_frac_first = Mole_flow_first / np.sum(Mole_flow_first)
        dpn = (np.sum(polymer_Mole_flow_first) / np.sum(polymer_Mole_flow_zeroth))

        c_idx_list = []
        for c_idx, comp in enumerate(GlobalComponentManager.component_list):
            if type(comp) is Component:
                c_idx_list.append(c_idx)

        mole_flow_properties = np.append(Mole_flow_zeroth[c_idx_list], pyo.value(np.sum(polymer_Mole_flow_zeroth)))
        mole_frac_properties = mole_flow_properties / np.sum(mole_flow_properties)

        pc_ftr_polymer = self.PropertiesMethod.param.r[-1]

        self.PropertiesMethod.param.m[-1] = pc_ftr_polymer * pyo.value(dpn) * self.PropertiesMethod.param.MW[-1]

        vm_liq = self.PropertiesMethod.calculate_molar_density_mixture(self.Temperature, self.Pressure,
                                                                       self.PropertiesMethod.param,
                                                                       mole_frac_properties, mole_frac_properties[-1],
                                                                       dpn)
        # molar density [mol/m3]

        q_out = np.sum(Mole_flow_first) / vm_liq

        concentration = np.array([Outflow.comp_dict[c]['mole_flow'] / q_out for c in Outflow.comp_dict])

        for f in self.Inflow.comp_dict:
            eq = self.Inflow.comp_dict[f]['mole_flow'] - Outflow.comp_dict[f][
                'mole_flow'] + self.ReactionSet.calculate_rate(f,
                                                               concentration) * self.Volume * 3600.
            mb_eq.append(eq)

        return mb_eq
