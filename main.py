import componentmanager
from component import *
from reactor import CstrSingleLiqPhase
from reactions import *
from flow import *
from proptiesmethod import *
import numpy as np

component_list = [Component("IB", '115-11-7', CompType.conventional), Component("IP", '78-79-5', CompType.conventional),
                  Component("HCL", '7647-01-0', CompType.conventional),
                  Component("EADC", '563-43-9', CompType.conventional),
                  Component("HEXANE", '110-54-3', CompType.conventional),
                  Component("CH3CL", '110-54-3', CompType.conventional)]

componentmanager.GlobalComponentManager.component_list_gen(component_list, 'CATION', 3)

print(len(componentmanager.GlobalComponentManager.component_list))
print(GlobalComponentManager.component_list)
reaction_set_1 = ReactionSet()

reaction_set_1.source_define(0, [0, 6], [1, 1], [1, 1], True,k_constant=3000)
reaction_set_1.source_define(0, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 10], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 11], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 12], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 10], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 11], [1, 1], [1, 1], True)
reaction_set_1.source_define(0, [0, 12], [1, 1], [1, 1], True)

reaction_set_1.source_define(2, [2, 3], [1, 1], [1, 1], True)
reaction_set_1.source_define(3, [2, 3], [1, 1], [1, 1], True)

reaction_set_1.source_define(6, [2, 3], [1, 1], [1, 1], False)
reaction_set_1.source_define(6, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(6, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(6, [0, 6], [1, 1], [1, 1], True)
reaction_set_1.source_define(6, [10, 15], [1, 1], [1, 1], False)
reaction_set_1.source_define(6, [11, 15], [1, 1], [1, 1], False)
reaction_set_1.source_define(6, [12, 15], [1, 1], [1, 1], False)

reaction_set_1.source_define(7, [0, 6], [1, 1], [1, 1], False)
reaction_set_1.source_define(7, [0, 7], [1, 1], [1, 1], True)

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


