from enum import Enum


class CompType(Enum):
    conventional = 0,
    polymer = 1,
    segment = 2


class Component:

    def __init__(self, formular, cas, comp_type, comp_mw=0.):
        self.name = formular
        self.CAS = cas
        self.type = comp_type
        self.polymer_mole_flow_relative = -2
        self.MW = comp_mw


class Polymer(Component):
    def __init__(self, segments, formular, cas, comp_type):
        super().__init__(formular, cas, comp_type)
        self.segment = segments
        self.type = CompType.polymer


class Segment(Component):

    def __init__(self, formular, cas, comp_type):
        super().__init__(formular, cas, comp_type)
