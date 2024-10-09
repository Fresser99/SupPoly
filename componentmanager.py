from component import *
import pandas as pd


class SpecificComponent:
    def __init__(self, name, polymer_flow_flag):
        self.polymer_mole_flow_relative = polymer_flow_flag
        self.name = name



class GlobalComponentManager:
    component_list = []

    @classmethod
    def component_list_gen(cls, clist: list[Component], poly_Type, poly_site):
        for c in clist:
            if c.type == CompType.conventional:
                cls.component_list.append(c)

        species_polymer_config = pd.read_csv('ReactConfig.csv')
        name_list = species_polymer_config.loc[:, poly_Type].values
        for idx, cname in enumerate(name_list):
            is_site_dependent = True if len(str(cname).split(';')) >= 2 and str(cname).split(';')[
                1] == 'NSite' else False

            is_mole_flow_relative = True if len(str(cname).split(';')) == 3 else False

            if is_site_dependent:
                for site in range(1, poly_site + 1):
                    name = str(cname).split(';')[0] + '[' + str(site) + ']'
                    if is_mole_flow_relative:
                        s_c = SpecificComponent(name, int(str(cname).split(';')[2]))
                    else:
                        s_c = SpecificComponent(name, -1)
                    cls.component_list.append(s_c)
            else:
                s_c = SpecificComponent(cname, -1)
                cls.component_list.append(s_c)
