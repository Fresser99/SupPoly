import numpy as np

import componentmanager
from flow import *
from componentmanager import GlobalComponentManager
from reactions import *
from proptiesmethod import *
import pyomo.environ as pyo


class CstrSingleLiqPhase:

    def __init__(self, t, p, v, inflow: Flow, rx_set: ReactionSet, prop: PropertiesMethod, q_spec=0):

        self.ReactionSet = rx_set
        self.Inflow = inflow
        self.PropertiesMethod = prop
        self.Temperature = t
        self.Pressure = p
        self.Volume = v
        self.q_spec = q_spec

    def mass_balance(self, Outflow: Flow):

        mb_eq = []

        polymer_Mole_flow_zeroth = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole][
                                                                    'polymer_flow_momentum'] == 0 else 0.0
             for mole in Outflow.comp_dict])

        polymer_Mole_flow_first = np.array(
            [pyo.value(Outflow.comp_dict[mole]['mole_flow']) if Outflow.comp_dict[mole][
                                                                    'polymer_flow_momentum'] == 1 else 0.0
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
        if self.q_spec != 0:
            q_out = np.sum(Mole_flow_first) * 1000 / self.q_spec
        else:
            q_out = np.sum(Mole_flow_first) * 1000 / vm_liq

        concentration = np.array([Outflow.comp_dict[c]['mole_flow'] / q_out for c in Outflow.comp_dict])
        scaling_factors = {}
        for f in self.Inflow.comp_dict:
            # 计算每个反应速率的预估值
            rate_estimate = abs(self.ReactionSet.calculate_rate(f, concentration)) * self.Volume * 3600
            # 计算缩放因子，这里使用最大值的倒数，防止除零
            inflow = abs(self.Inflow.comp_dict[f]['mole_flow'])
            outflow = abs(Outflow.comp_dict[f]['mole_flow'])
            max_value = max(pyo.value(inflow), pyo.value(outflow), pyo.value(rate_estimate), 1e-8)
            scaling_factors[f] = 1.0 / max_value
        for f in self.Inflow.comp_dict:
            # 原始方程
            eq = (self.Inflow.comp_dict[f]['mole_flow'] - Outflow.comp_dict[f]['mole_flow'] +
                  self.ReactionSet.calculate_rate(f, concentration) * self.Volume * 3600)
            # 应用缩放因子
            scaled_eq = scaling_factors[f] * eq
            mb_eq.append(scaled_eq)

        return mb_eq
        mb_eq.append(eq)

        return mb_eq


class PFRSingleliqPhase:

    def __init__(self, t, p, l, D, inflow: Flow, rx_set: ReactionSet, prop: PropertiesMethod):
        self.ReactionSet = rx_set
        self.Inflow = inflow
        self.PropertiesMethod = prop
        self.Temperature = t
        self.Pressure = p
        self.Area = D
        self.Length = l

    def compute_dpn(self, model, z):
        polymer_zeroth_mole_flow = sum(
            model.F[c, z] if self.Inflow.comp_dict[c]['polymer_flow_momentum'] == 0 else 0.0 for c in
            self.Inflow.comp_dict)

        polymer_first_mole_flow = sum(
            model.F[c, z] if self.Inflow.comp_dict[c]['polymer_flow_momentum'] == 1 else 0.0 for c in
            self.Inflow.comp_dict)

        mole_flow_zeroth = sum(
            model.F[c, z] if self.Inflow.comp_dict[c]['polymer_flow_momentum'] in [0, -2] else 0.0 for c in
            self.Inflow.comp_dict)

        mole_flow_first = sum(
            model.F[c, z] if self.Inflow.comp_dict[c]['polymer_flow_momentum'] in [1, -2] else 0.0 for c in
            self.Inflow.comp_dict)

        return polymer_first_mole_flow / polymer_zeroth_mole_flow

    def volume_flow_rate_rule(self, model, z):

        polymer_zeroth_mole_flow = np.array(
            [pyo.value(model.F[c, z]) if self.Inflow.comp_dict[c]['polymer_flow_momentum'] == 0 else 0.0
             for c in
             self.Inflow.comp_dict])

        mole_flow_zeroth = np.array([pyo.value(model.F[c, z]) if self.Inflow.comp_dict[c][
                                                                     'polymer_flow_momentum'] in [0,
                                                                                                  -2] else 0.0
                                     for c in
                                     self.Inflow.comp_dict])

        mole_flow_first = [
            pyo.value(model.F[c, z]) if self.Inflow.comp_dict[c]['polymer_flow_momentum'] in [1,
                                                                                              -2] else 0.0
            for c in
            self.Inflow.comp_dict]

        c_idx_list = []
        for c_idx, comp in enumerate(GlobalComponentManager.component_list):
            if type(comp) is Component:
                c_idx_list.append(c_idx)

        mole_flow_properties = np.append(mole_flow_zeroth[c_idx_list], np.sum(polymer_zeroth_mole_flow))
        mole_frac_properties = mole_flow_properties / np.sum(mole_flow_properties)

        pc_ftr_polymer = self.PropertiesMethod.param.r[-1]

        self.PropertiesMethod.param.m[-1] = pc_ftr_polymer * pyo.value(model.dpn[z]) * \
                                            self.PropertiesMethod.param.MW[-1]

        vm_liq = self.PropertiesMethod.calculate_molar_density_mixture(self.Temperature,
                                                                       self.Pressure,
                                                                       self.PropertiesMethod.param,
                                                                       mole_frac_properties,
                                                                       mole_frac_properties[-1],
                                                                       model.dpn[z]) / 1000

        return vm_liq

    def concentration_rule(self, model, comp, z):
        return model.F[comp, z] / model.V_flow[z]

    def mass_balance(self, model, comp, z):

        concentrations = [model.F[c.name, z] / model.V_flow[z] for c in GlobalComponentManager.component_list]
        rate = self.ReactionSet.calculate_rate(comp, concentrations)
        return model.dFdz[comp, z] == rate * self.Area * 3600
