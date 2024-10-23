from flow import *
import pyomo.environ as pyo


class Mixer:
    def __init__(self, inlet_flow_list):

        self.split_Out_flow = Flow(t=101, p=101, name='mix')
        self.inlet_flow_list = inlet_flow_list

    def mass_balance(self):

        for c in self.split_Out_flow.comp_dict:
            mix_flow = 0.0
            for idx_flow, flow in enumerate(self.inlet_flow_list):
                mix_flow = mix_flow + flow.comp_dict[c]['mole_flow']
            self.split_Out_flow.comp_dict[c]['mole_flow'] = mix_flow


class Spliter:

    def __init__(self, inlet_flow: Flow, split_frac):

        self.split_Out_flow_list = [Flow(t=101, p=0, name='split_flow_out_' + str(i)) for i in range(len(split_frac))]
        self.split_frac = split_frac
        self.inlet_flow = inlet_flow

    def mass_balance(self):

        for flow_idx, flow in enumerate(self.split_Out_flow_list):
            for comp in self.inlet_flow.comp_dict:
                flow.comp_dict[comp]['mole_flow'] = self.inlet_flow.comp_dict[comp]['mole_flow']* \
                                                    self.split_frac[
                                                        flow_idx]
