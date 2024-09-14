from component import Component


class GlobalComponentManager:

    component_list=[]
    @classmethod
    def component_list_gen(cls, clist: list[Component]):

        for c in clist:
            cls.component_list.append(c.Formular)

        cls.component_list.append('ion-pair')
        cls.component_list.append('active_species')
        cls.component_list.append('p1_ion')
        cls.component_list.append('zeroth_mom_live')
        cls.component_list.append('first_mom_live')
        cls.component_list.append('second_mom_live')
        cls.component_list.append('zeroth_mom_dead')
        cls.component_list.append('first_mom_dead')
        cls.component_list.append('second_mom_dead')
        cls.component_list.append('counter_ion')

