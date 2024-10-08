from component import Component
import pandas as pd


class GlobalComponentManager:
    component_list = []

    @classmethod
    def component_list_gen(cls, clist: list[Component],poly_Type,poly_site):
        for c in clist:
            cls.component_list.append(c.Formular)

        species_polymer_config = pd.read_csv('ReactConfig.csv')
        name_list = species_polymer_config.loc[:, poly_Type].values
        for idx, cname in enumerate(name_list):
            is_site_dependent = True if len(str(cname).split(';')) == 2 and str(cname).split(';')[
                1] == 'NSite' else False

            if is_site_dependent:
                for site in range(1,poly_site+1):
                    cls.component_list.append(str(cname).split(';')[0]+'['+str(site)+']')
            else:
                cls.component_list.append(cname)
