from componentmanager import GlobalComponentManager
from component import Component
import numpy as np
from enum import Enum


class MetaComponentType(Enum):
    CATALYST = 1,
    INITIATOR = 2,
    MONOMER = 3,


class MetaReactionType(Enum):
    ACTIVATE = 1,
    INITIATION = 2,
    PROPAGATION = 3,
    TRANSFER = 4,
    TERMINATION = 5,
    TOXIFICATION = 6


class ReactionSet:

    def __init__(self):
        self.source_dict = {c.name: [] for c in GlobalComponentManager.component_list}

    def output_rx_matrix(self):

        for r_idx, r in enumerate(self.reactions_set):

            for c_idx, c in enumerate(GlobalComponentManager.component_list):

                for cr in r.lf:

                    if c == cr:
                        self.Rx_matrix[r_idx, c_idx] = -1
                for cr in r.rt:

                    if c == cr:
                        self.Rx_matrix[r_idx, c_idx] = 1

    def print_kinetics(self):

        rate = {c: [] for c in GlobalComponentManager.component_list}
        for c_idx, c in enumerate(GlobalComponentManager.component_list):

            for r_idx, r in enumerate(self.reactions_set):

                if self.Rx_matrix[r_idx, c_idx] == -1.:

                    reactant = [index for index, value in enumerate(self.Rx_matrix[r_idx, :]) if value == -1]

                    if len(reactant) == 2:
                        rate[c].append(
                            '-' + '[' + GlobalComponentManager.component_list[reactant[0]] + ']' + '*' + '[' +
                            GlobalComponentManager.component_list[reactant[1]] + ']')
        return rate

    def set_comp_type(self, comp: Component, type: MetaComponentType):

        if type == MetaComponentType.CATALYST:

            self.configure_dict["CAT"].append(comp)
        elif type == MetaComponentType.INITIATOR:
            self.configure_dict["INITIATOR"].append(comp)
        elif type == MetaComponentType.MONOMER:
            self.configure_dict["MONOMER"].append(comp)

    def source_define(self, idx_obj, powerlaw_idx, k_constant, param=None, order=None, is_sink=True, ):

        if order is None:
            order = [1, 1]
        if param is None:
            param = [1, 1]
        cc = GlobalComponentManager.component_list[idx_obj]
        comp = cc.name
        # comp = GlobalComponentManager.component_list[idx_obj].Formular if GlobalComponentManager.component_list[idx_obj] is Component else GlobalComponentManager.component_list[idx_obj].name

        source = [i for i in powerlaw_idx]
        source = source + param + order
        if is_sink:
            source.append(-1)
        else:

            source.append(1)
        source.append(k_constant['value'])
        source.append(k_constant['name'])
        self.source_dict[comp].append(source)

    def calculate_rate(self, comp_key, concen):

        tot_rate = 0.
        for c in self.source_dict[comp_key]:
            tot_rate = tot_rate + c[7] * concen[c[0]] * c[2] * concen[c[1]] * c[3] * c[6]

        return tot_rate

    def preview_reaction_equations(self):

        for i in self.source_dict:
            eq = 'r' + '[' + str(i) + ']' + '='
            for s in self.source_dict[i]:
                eq = eq + ('-' if s[6] == -1 else '+') + str(s[8]) + '*' + str(
                    s[2]) + '*' + '[' + GlobalComponentManager.component_list[s[0]].name + ']' + '*' + str(
                    s[3]) + '*' + '[' + GlobalComponentManager.component_list[s[1]].name + ']'
            print(eq)


class MetaReaction:

    def __init__(self, left_side_comps, right_side_comps, is_reservable):
        self.lf = left_side_comps
        self.rt = right_side_comps
        self.is_reserve = is_reservable
