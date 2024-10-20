import numpy as np

from componentmanager import *


class Flow:

    def __init__(self, t, p, name):
        self.Temperature = t
        self.Pressure = p
        self.Name = name
        self.comp_dict = {
            c.name: {"mass_flow": 0., "mole_flow": 0., "polymer_flow_momentum": c.polymer_mole_flow_relative} for c in
            GlobalComponentManager.component_list}

    def get_mole_frac(self):
        tot_mole = []
        for c in self.comp_dict:
            tot_mole.append(self.comp_dict[c]['mole_flow'])

        tot = np.sum(tot_mole)
        print(tot)

    def set_MassFlow_conventional(self, flow_arr):
        for c in flow_arr:
            assert c in self.comp_dict, '组分：' + '\'' + str(c) + '\'' + '不在组分列表里'

            self.comp_dict[c]['mass_flow'] = flow_arr[c]
        self.mass_2_mole()

    def mass_2_mole(self):

        for idx, c in enumerate(self.comp_dict):
            self.comp_dict[c]['mole_flow'] = self.comp_dict[c]['mass_flow'] / GlobalComponentManager.component_list[
                idx].MW if self.comp_dict[c]['polymer_flow_momentum'] == -2 else 0.0

