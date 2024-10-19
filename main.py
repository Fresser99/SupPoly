from pyomo.opt import SolverFactory

import componentmanager
from component import *
from reactor import CstrSingleLiqPhase
from reactions import *
from flow import *
from proptiesmethod import *
import numpy as np
from pcsaft.param import *
from pcsaft.pcsaft import *
import pyomo.environ as pyo

site_num = 3

IIR = Polymer([Segment('IB-R', 'IB-R-seg', CompType.segment)], 'IIR', 'IIR', CompType.polymer)

component_list = [Component("IB", '115-11-7', CompType.conventional, 56.1075),
                  Component("IP", '78-79-5', CompType.conventional, 68.1185),
                  Component("HCL", '7647-01-0', CompType.conventional, 36.4606),
                  Component("EADC", '563-43-9', CompType.conventional, 126.949),
                  Component("HEXANE", '110-54-3', CompType.conventional, 86.1772),
                  Component("CH3CL", '74-87-3', CompType.conventional, 50.4875), IIR]

componentmanager.GlobalComponentManager.component_list_gen(component_list, 'CATION', site_num)

param = Param
properties_package = IIR_PCSAFT(param)
param.m = np.array([])
param.e = np.array([])
param.s = np.array([])
param.r = np.array([])
param.MW = np.array([])

for c in GlobalComponentManager.component_list:

    if type(c) is Component:
        param.m = np.append(param.m, np.float32(properties_package.pcsaft.retrive_param_from_DB(c.CAS, 'PCFTM')))
        param.e = np.append(param.e, np.float32(properties_package.pcsaft.retrive_param_from_DB(c.CAS, 'PCFTU')))
        param.s = np.append(param.s, np.float32(properties_package.pcsaft.retrive_param_from_DB(c.CAS, 'PCFTV')))
        param.MW = np.append(param.MW, np.float32(properties_package.pcsaft.retrive_param_from_DB(c.CAS, 'MW')))

param.m = np.append(param.m, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IP-R', 'PCFTR')))
param.e = np.append(param.e, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IP-R', 'PCFTU')))
param.s = np.append(param.s, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IP-R', 'PCFTV')))
param.MW = np.append(param.MW, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IP-R', 'MW')))
param.r = np.append(param.r, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IP-R', 'PCFTR')))

param.k_ij = np.zeros([len(component_list), len(component_list)])
# print(param.m)
# print(param.s)
# print(param.e)
# print(properties_package.param.m)
# print(len(componentmanager.GlobalComponentManager.component_list))
# for c in GlobalComponentManager.component_list:
#     print(c.name)
reaction_set_1 = ReactionSet()
kinetic = {'ka': 0.0019, 'ki(1)': 0.001, 'ki(2)': 0, 'ki(3)': 0.0, 'kp(1)': 996, 'kp(2)': 0, 'kp(3)': 0.0,
           'ktm(1)': 0.28, 'ktm(2)': 0.0, 'ktm(3)': 0.0, 'kd(1)': 0.04, 'kd(2)': 0.0, 'kd(3)': 0.0}

# IB + ion-pair ——> P1[1] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[2] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[3] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1] ——>  Y0[1]
reaction_set_1.source_define(0, [0, 10], {'name': 'Kp(1)', 'value': kinetic['kp(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2] ——>  Y0[2]
reaction_set_1.source_define(0, [0, 11], {'name': 'Kp(2)', 'value': kinetic['kp(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3] ——>  Y0[3]
reaction_set_1.source_define(0, [0, 12], {'name': 'Kp(3)', 'value': kinetic['kp(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1] ——> D0[1]+P1[1]
reaction_set_1.source_define(0, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2] ——> D0[2]+P1[1]
reaction_set_1.source_define(0, [0, 11], {'name': 'Ktm(1)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3] ——> D0[3]+P1[3]
reaction_set_1.source_define(0, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# HCL + EADC ——> ion-pair
reaction_set_1.source_define(2, [2, 3], {'name': 'Ka', 'value': kinetic['ka']}, [1, 1], [1, 1], True)
# EADC + HCL ——> ion-pair
reaction_set_1.source_define(3, [2, 3], {'name': 'Ka', 'value': kinetic['ka']}, [1, 1], [1, 1], True)

# HCL + EADC ——> ion-pair
reaction_set_1.source_define(6, [2, 3], {'name': 'Ka', 'value': kinetic['ka']}, [1, 1], [1, 1], False)

# IB + ion-pair ——> P1[1] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], True)

# IB + ion-pair ——> P1[2] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], True)

# IB + ion-pair ——> P1[3] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], True)

# Y0[1] + conuter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(6, [10, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], False)
# Y0[2] + conuter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(6, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + conuter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(6, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# IB + ion-pair ——> p1_ion[1] + counter-ion
reaction_set_1.source_define(7, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], False)
# IB + ion-pair ——> p1_ion[2] + counter-ion
reaction_set_1.source_define(8, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], False)
# IB + ion-pair ——> p1_ion[3] + counter-ion
reaction_set_1.source_define(9, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], False)

# IB + p1_ion[1] ——> p2[1]
reaction_set_1.source_define(7, [0, 7], {'name': 'Kp(1)', 'value': kinetic['kp(1)']}, [1, 1], [1, 1], True)
# IB + p1_ion[2] ——> p2[2]
reaction_set_1.source_define(8, [0, 8], {'name': 'Kp(2)', 'value': kinetic['kp(2)']}, [1, 1], [1, 1], True)
# IB + p1_ion[3] ——> p2[3]
reaction_set_1.source_define(9, [0, 9], {'name': 'Kp(3)', 'value': kinetic['kp(3)']}, [1, 1], [1, 1], True)

# p1_ion[1] + IB ——> D1[1]: p1 monomer transfer
reaction_set_1.source_define(7, [0, 7], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> D1[2]
reaction_set_1.source_define(8, [0, 8], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> D1[3]
reaction_set_1.source_define(9, [0, 9], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# p1_ion[1] + counter_ion ——> D1[1] + ion-pair
reaction_set_1.source_define(7, [7, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], True)
# p1_ion[2] + counter_ion ——> D1[2] + ion-pair
reaction_set_1.source_define(8, [8, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + counter_ion ——> D1[3] + ion-pair
reaction_set_1.source_define(9, [9, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(7, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2]
reaction_set_1.source_define(8, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(9, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# ion-pair + IB ——> P1_ion[1] ( Y0[1] )
reaction_set_1.source_define(10, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> P1_ion[2] ( Y0[2] )
reaction_set_1.source_define(11, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> P1_ion[3] ( Y0[3] )
reaction_set_1.source_define(12, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(10, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + pi_ion[2]
reaction_set_1.source_define(11, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(12, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(10, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# Y0[2] + IB ——> D0[2] + pi_ion[2]
reaction_set_1.source_define(11, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# Y0[3] + IB ——> D0[3] + pi_ion[3]
reaction_set_1.source_define(12, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(10, [10, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], True)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(11, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(12, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] (1*Y0[1]=Y1[1] n=1)
reaction_set_1.source_define(13, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] (1*Y0[2]=Y1[2] n=1)
reaction_set_1.source_define(14, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] (1*Y0[3]=Y1[3] n=1)
reaction_set_1.source_define(15, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1](1*Y0[1]=Y1[1] n=1)
reaction_set_1.source_define(13, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2](1*Y0[2]=Y1[2] n=1)
reaction_set_1.source_define(14, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3](1*Y0[3]=Y1[3] n=1)
reaction_set_1.source_define(15, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y1[1] + IB ——> Y1[1]
reaction_set_1.source_define(13, [0, 10], {'name': 'Kp(1)', 'value': kinetic['kp(1)']}, [1, 1], [1, 1], False)
# Y1[2] + IB ——> Y1[2]
reaction_set_1.source_define(14, [0, 11], {'name': 'Kp(2)', 'value': kinetic['kp(2)']}, [1, 1], [1, 1], False)
# Y1[3] + IB ——> Y1[3]
reaction_set_1.source_define(15, [0, 12], {'name': 'Kp(3)', 'value': kinetic['kp(3)']}, [1, 1], [1, 1], False)

# Y1[1] + IB ——> Dn[1] + p1_ion[1]
reaction_set_1.source_define(13, [0, 13], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# Y1[2] + IB ——> Dn[2] + p1_ion[2]
reaction_set_1.source_define(14, [0, 14], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# Y1[3] + IB ——> Dn[3] + p1_ion[3]
reaction_set_1.source_define(15, [0, 15], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# Y1[1] + counter-ion ——> Dn[1] + ion-pair[1]
reaction_set_1.source_define(13, [13, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], True)
# Y1[2] + counter-ion ——> Dn[2] + ion-pair[2]
reaction_set_1.source_define(14, [14, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y1[3] + counter-ion ——> Dn[3] + ion-pair[3]
reaction_set_1.source_define(15, [15, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] (1^2*Y0[1]=Y2[1] n=1)
reaction_set_1.source_define(16, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] (1^2*Y0[2]=Y2[2] n=1)
reaction_set_1.source_define(17, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] (1^2*Y0[3]=Y2[3] n=1)
reaction_set_1.source_define(18, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1](1^2*Y0[1]=Y2[1] n=1)
reaction_set_1.source_define(16, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2](1^2*Y0[1]=Y2[2] n=1)
reaction_set_1.source_define(17, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3](1^2*Y0[1]=Y2[3] n=1)
reaction_set_1.source_define(18, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# 2Y1[1] + IB ——> Y2[1]
reaction_set_1.source_define(16, [0, 13], {'name': 'Kp(1)', 'value': kinetic['kp(1)']}, [1, 2], [1, 1], False)
# 2Y1[2] + IB ——> Y2[2]
reaction_set_1.source_define(17, [0, 14], {'name': 'Kp(2)', 'value': kinetic['kp(2)']}, [1, 2], [1, 1], False)
# 2Y1[3] + IB ——> Y2[3]
reaction_set_1.source_define(18, [0, 15], {'name': 'Kp(3)', 'value': kinetic['kp(3)']}, [1, 2], [1, 1], False)

# Y0[1] + IB ——> Y2[1]
reaction_set_1.source_define(16, [0, 10], {'name': 'Kp(1)', 'value': kinetic['kp(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> Y2[2]
reaction_set_1.source_define(17, [0, 11], {'name': 'Kp(2)', 'value': kinetic['kp(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> Y2[3]
reaction_set_1.source_define(18, [0, 12], {'name': 'Kp(3)', 'value': kinetic['kp(3)']}, [1, 1], [1, 1], False)

# Y2[1] + IB ——> D2[1] + p1_ion[1]
reaction_set_1.source_define(16, [0, 16], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# Y2[2] + IB ——> D2[1] + p1_ion[2]
reaction_set_1.source_define(17, [0, 17], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# Y2[3] + IB ——> D2[3] + p1_ion[3]
reaction_set_1.source_define(18, [0, 18], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# Y2[1] +counter-ion+ ——> D2[1] + ion-pair[1]
reaction_set_1.source_define(16, [16, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], True)
# Y2[2] +counter-ion+ ——> D2[2] + ion-pair[2]
reaction_set_1.source_define(17, [17, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y2[3] +counter-ion+ ——> D2[3] + ion-pair[3]
reaction_set_1.source_define(18, [18, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(19, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2]
reaction_set_1.source_define(20, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(21, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(19, [10, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], False)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(20, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(21, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(19, [0, 7], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(20, [0, 8], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(21, [0, 9], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(22, [0, 13], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(23, [0, 14], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(24, [0, 15], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y1[1] + counter_ion ——> D1[1] +ion-pair[1]
reaction_set_1.source_define(22, [13, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], False)
# Y1[2] + counter_ion ——> D1[2] +ion-pair[2]
reaction_set_1.source_define(23, [14, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[3] + counter_ion ——> D1[3] +ion-pair[3]
reaction_set_1.source_define(24, [15, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(22, [0, 7], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(23, [0, 8], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(24, [0, 9], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# Y2[1] + IB ——> D2[1] + p1_ion[1]
reaction_set_1.source_define(25, [0, 22], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y2[2] + IB ——> D2[2] + p1_ion[2]
reaction_set_1.source_define(26, [0, 23], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y2[3] + IB ——> D2[3] + p1_ion[3]
reaction_set_1.source_define(27, [0, 24], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y2[1] + counter_ion ——> D2[1] +ion-pair[1]
reaction_set_1.source_define(25, [16, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], False)
# Y2[2] + counter_ion ——> D2[2] +ion-pair[2]
reaction_set_1.source_define(26, [17, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y2[3] + counter_ion ——> D2[3] +ion-pair[3]
reaction_set_1.source_define(27, [18, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(25, [0, 7], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(26, [0, 8], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(27, [0, 9], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(1)', 'value': kinetic['ki(1)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(2)', 'value': kinetic['ki(2)']}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(3)', 'value': kinetic['ki(3)']}, [1, 1], [1, 1], False)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(28, [10, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], True)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(28, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(28, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

reaction_set_1.preview_reaction_equations()
flow_toR130 = Flow(100, 103, 'to_R130')
flow_toR130.set_MassFlow_conventional(
    {'IB': 4070., 'IP': 0, 'HCL': 0.47, 'EADC': 12, 'HEXANE': 150, 'CH3CL': 8600})
flow_outR130 = Flow(t=0, p=101, name='R130_out')
reactor = CstrSingleLiqPhase(198., 10132500., 24, flow_toR130, reaction_set_1, properties_package)

model = pyo.ConcreteModel()

comp_names = [c.name for c in GlobalComponentManager.component_list]
init_values = {v: flow_toR130.comp_dict[v]['mole_flow'] + 1e-6 for v in flow_toR130.comp_dict}
model.outflows = pyo.Var(comp_names, initialize=init_values, domain=pyo.NonNegativeReals)

for c in flow_outR130.comp_dict:
    flow_outR130.comp_dict[c]['mole_flow'] = model.outflows[c]

model.eqs = pyo.ConstraintList()

for eq in reactor.mass_balance(flow_outR130):
    model.eqs.add(eq == 0.)

solver = SolverFactory('ipopt')
results = solver.solve(model, tee=True)

if (results.solver.status == pyo.SolverStatus.ok) and (
        results.solver.termination_condition == pyo.TerminationCondition.optimal):
    print("Optimal solution found.")
    for name in comp_names:
        print(f"{name}: {pyo.value(model.outflows[name])}")
    print("MWN: " + str(
        (pyo.value(model.outflows['first_mom_live[1]'])+pyo.value(model.outflows['first_mom_dead[1]'])) / (pyo.value(model.outflows['zeroth_mom_dead[1]'])+pyo.value(model.outflows['zeroth_mom_dead[1]'])) * 56))
    print("IB进料量：" + str(flow_toR130.comp_dict['IB']['mole_flow']))
    print("转化率: " + str(
        (flow_toR130.comp_dict['IB']['mole_flow'] - pyo.value(model.outflows['IB'])) / flow_toR130.comp_dict['IB'][
            'mole_flow'] * 100))
else:
    print("Solver did not find an optimal solution.")

print(str(pyo.value(model.outflows['first_mom_live[1]']) + pyo.value(model.outflows['first_mom_dead[1]'])))
