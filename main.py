import componentmanager
from component import Component
from reactor import CstrSingleLiqPhase
from reactions import *
from flow import *
from proptiesmethod import *

component_list = [Component("IB"), Component("IP"), Component("HCL"), Component("EADC"), Component("HEXANE"),
                  Component("CH3CL")]

componentmanager.GlobalComponentManager.component_list_gen(component_list)

# Rx_asso=Meta_Reaction(['HCL','EADC'],['ion-pair'],False)
# Rx_active=Meta_Reaction(['ion-pair','IB'],['active_species','counter_ion'],False)
# Rx_init=Meta_Reaction(['active_species','IB'],['p1_ion'],False)
# Rx_prop=Meta_Reaction(['zeroth_mom_live','IB'],['zeroth_mom_live'],False)
# Rx_chat=Meta_Reaction(['zeroth_mom_live','IB'],['zeroth_mom_dead','active_species'],False)
# Rx_term=Meta_Reaction(['zeroth_mom_live','counter_ion'],['zeroth_mom_dead','ion-pair'],False)
# # Rx_term=Meta_Reaction(['first_mom_live','IB'],[''])
reaction_set_1 = ReactionSet()
# print(reaction_set_1.Rx_matrix)
# ss=reaction_set_1.pirnt_kinetics()
# print(ss)
print(GlobalComponentManager.component_list)
reaction_set_1.source_define(0, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 7], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 9], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 9], [1, 1], [1, 1], True)

reaction_set_1.source_define(2, [2, 3], [1, 1], [1, 1], True)
reaction_set_1.source_define(3, [2, 3], [1, 1], [1, 1], True)

reaction_set_1.source_define(6, [2, 3], [1, 1], [1, 1], False)
reaction_set_1.source_define(6, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(6, [9, 15], [1, 1], [1, 1], False)

reaction_set_1.source_define(7, [0, 6], [1, 1], [1, 1], False)
reaction_set_1.source_define(7, [0, 9], [1, 1], [1, 1], False)

reaction_set_1.source_define(8, [0, 7], [1, 1], [1, 1], False)
reaction_set_1.source_define(8, [0, 8], [1, 1], [1, 1], True)

reaction_set_1.source_define(9, [0, 7], [1, 1], [1, 1], False)
reaction_set_1.source_define(9, [0, 9], [1, 1], [1, 1], True)
reaction_set_1.source_define(9, [9, 15], [1, 1], [1, 1], True)

reaction_set_1.source_define(10, [0, 7], [1, 1], [1, 1], False)
reaction_set_1.source_define(10, [0, 9], [1, 1], [1, 1], False)
reaction_set_1.source_define(10, [0, 9], [1, 1], [1, 1], True)
reaction_set_1.source_define(10, [9, 15], [1, 1], [1, 1], True)

reaction_set_1.source_define(11, [0, 7], [1, 1], [1, 1], False)
reaction_set_1.source_define(11, [0, 10], [1, 2], [1, 1], False)
reaction_set_1.source_define(11, [0, 9], [1, 1], [1, 1], False)
reaction_set_1.source_define(11, [0, 11], [1, 1], [1, 1], True)
reaction_set_1.source_define(11, [15, 11], [1, 1], [1, 1], True)

reaction_set_1.source_define(12, [0, 9], [1, 1], [1, 1], False)
reaction_set_1.source_define(12, [15, 9], [1, 1], [1, 1], False)

reaction_set_1.source_define(13, [0, 10], [1, 1], [1, 1], False)
reaction_set_1.source_define(13, [15, 10], [1, 1], [1, 1], False)

reaction_set_1.source_define(14, [0, 11], [1, 1], [1, 1], False)
reaction_set_1.source_define(14, [15, 11], [1, 1], [1, 1], False)

reaction_set_1.source_define(15, [0, 6], [1, 1], [1, 1], False)
reaction_set_1.source_define(15, [9, 15], [1, 1], [1, 1], True)

flow_toR130 = Flow(100, 103, 'to_R130')

properties_method = UserMethod()
reactor = CstrSingleLiqPhase(100., 100., 30, flow_toR130, reaction_set_1, properties_method)

print(reaction_set_1.source_dict)

print(flow_toR130.comp_dict)
flow_toR130.get_mole_frac()
T = 178.15  # K
P = 1e6  # Pa
Tc = np.array([190.6, 305.4, 369.8, 425.2])  # K
Pc = np.array([4.6e6, 4.88e6, 4.25e6, 3.8e6])  # Pa
omega = np.array([0.011, 0.099, 0.153, 0.199])
param_list={'Tc':Tc,'Pc':Pc,'Omega':omega}
y = np.array([0.3, 0.4, 0.2, 0.1])  # 摩尔分数
den,z=properties_method.calculate_molar_density_mixture(y,T,P,param_list)


print(den)

