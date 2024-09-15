import numpy as np

from componentmanager import *


class Flow:

    def __init__(self, t, p, name):
        self.Temperature = t
        self.Pressure = p
        self.Name = name
        self.comp_dict = {c: {"mass_flow": 0., "mole_flow": 0.} for c in GlobalComponentManager.component_list}

    def get_mole_frac(self):
        tot_mole = []
        for c in self.comp_dict:
            tot_mole.append(self.comp_dict[c]['mole_flow'])

        tot = np.sum(tot_mole)
        print(tot)
