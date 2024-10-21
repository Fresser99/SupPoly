from pyomo import dae
from pyomo.opt import SolverFactory

import componentmanager
from component import *
from reactor import *
from reactions import *
from flow import *
from proptiesmethod import *
import numpy as np
from pcsaft.param import *
from pcsaft.pcsaft import *
import pyomo.environ as pyo
from solvermanage import *
import matplotlib.pyplot as plt
from OperationUnit import *

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

param.m = np.append(param.m, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IB-R', 'PCFTR')))
param.e = np.append(param.e, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IB-R', 'PCFTU')))
param.s = np.append(param.s, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IB-R', 'PCFTV')))
param.MW = np.append(param.MW, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IB-R', 'MW')))
param.r = np.append(param.r, np.float32(properties_package.pcsaft.retrive_param_from_DB('seg-IB-R', 'PCFTR')))

param.k_ij = np.zeros([len(component_list), len(component_list)])
# print(param.m)
# print(param.s)
# print(param.e)
# print(properties_package.param.m)
# print(len(componentmanager.GlobalComponentManager.component_list))
# for c in GlobalComponentManager.component_list:
#     print(c.name)
reaction_set_1 = ReactionSet()
kinetic = {'ka': 0.0018, 'ki(1)': 0.0002, 'ki(2)': 0, 'ki(3)': 0.0, 'kp(1)': 981.1, 'kp(2)': 0, 'kp(3)': 0.0,
           'ktm(1)': 0.3, 'ktm(2)': 0.0, 'ktm(3)': 0.0, 'kd(1)': 0.01, 'kd(2)': 0.0, 'kd(3)': 0.0}

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

# reaction_set_1.preview_reaction_equations()

tear_flow = {cmp.name: 0.0 for cmp in GlobalComponentManager.component_list}
returnflow = {cmp.name: 0.0 for cmp in GlobalComponentManager.component_list}
flow_toR130 = Flow(100, 103, 'to_R130')
flow_toR130.set_MassFlow_conventional(
    {'IB': 4070., 'IP': 0, 'HCL': 0.47, 'EADC': 12, 'HEXANE': 150, 'CH3CL': 8600})



for comps in flow_toR130.comp_dict:
    flow_toR130.comp_dict[comps]['mole_flow'] += tear_flow[comps]

reactor = CstrSingleLiqPhase(198., 10132500., 3, flow_toR130, reaction_set_1, properties_package,q_spec=6000)
model1 = ReactorModel("R130", reactor, flow_toR130)

reactor2 = CstrSingleLiqPhase(198., 10132500., 3, model1.outlet_flow, reaction_set_1, properties_package)
model2 = ReactorModel("R130-2", reactor2, model1.outlet_flow)

reactor3 = CstrSingleLiqPhase(198., 10132500., 3, model2.outlet_flow, reaction_set_1, properties_package)
model3 = ReactorModel("R130-3", reactor3, model2.outlet_flow)

reactor4 = CstrSingleLiqPhase(198., 10132500., 3, model3.outlet_flow, reaction_set_1, properties_package)
model4 = ReactorModel("R130-4", reactor4, model3.outlet_flow)

solve_manager = SolverManager()
solve_manager.add_model(model1)
solve_manager.add_model(model2)
solve_manager.add_model(model3)
solve_manager.add_model(model4)
# solve_manager.add_model(tube_model)
result = solve_manager.solve_sequence()
PostProcess.process_results(result)
# IB_flow = []
# act_sp = []
# for i in result:
#     IB_flow.append(i['outflows']['IB'])
#     act_sp.append(i['outflows']['p1_ion[1]'])
# print(IB_flow)
# plt.plot(act_sp)
# plt.show()

# print(pyo.value(model4.outlet_flow.comp_dict['IB']['mole_flow']))

spliter_0 = Spliter(model4.outlet_flow, [0.001453043, 0.418959329, 0.274750957, 0.143179036, 0.161657635])
spliter_0.mass_balance()
tubeR = PFRSingleliqPhase(198, 10132500, 3.1, 0.18372636, spliter_0.split_Out_flow_list[1], reaction_set_1,
                          properties_package)
tube_model = ReactorModel_PFR('tube_outer_1', tubeR, inlet_flow=spliter_0.split_Out_flow_list[1])
solve_manager2 = SolverManager()
solve_manager2.add_model(tube_model)
res = solve_manager2.solve_sequence()
PostProcess.process_results(res)

tubeR_2 = PFRSingleliqPhase(198, 10132500, 3.1, 0.1836915, spliter_0.split_Out_flow_list[2], reaction_set_1,
                            properties_package)
tube_model_2 = ReactorModel_PFR('tube_outer_2', tubeR_2, inlet_flow=spliter_0.split_Out_flow_list[2])
discretizer = pyo.TransformationFactory('dae.finite_difference')
discretizer.apply_to(tube_model_2.model, nfe=50, scheme='BACKWARD')
solver2 = SolverFactory('ipopt')
solver2.solve(tube_model_2.model)
tubeR_3 = PFRSingleliqPhase(198, 10132500, 3.1, 0.12245611, spliter_0.split_Out_flow_list[3], reaction_set_1,
                            properties_package)
tube_model_3 = ReactorModel_PFR('tube_inner_1', tubeR_3, inlet_flow=spliter_0.split_Out_flow_list[3])
discretizer = pyo.TransformationFactory('dae.finite_difference')
discretizer.apply_to(tube_model_3.model, nfe=50, scheme='BACKWARD')
solver2 = SolverFactory('ipopt')
solver2.solve(tube_model_3.model)

tubeR_4 = PFRSingleliqPhase(198, 10132500, 3.1, 0.12246347, spliter_0.split_Out_flow_list[4], reaction_set_1,
                            properties_package)
tube_model_4 = ReactorModel_PFR('tube_inner_2', tubeR_4, inlet_flow=spliter_0.split_Out_flow_list[4])
discretizer = pyo.TransformationFactory('dae.finite_difference')
discretizer.apply_to(tube_model_4.model, nfe=50, scheme='BACKWARD')
solver2 = SolverFactory('ipopt')
solver2.solve(tube_model_4.model)

# plt.plot(pyo.value(tube_model.model.F['IB', :]))
# plt.plot(pyo.value(tube_model_2.model.F['IB', :]))
# plt.plot(pyo.value(tube_model_3.model.F['IB', :]))
# plt.plot(pyo.value(tube_model_4.model.F['IB', :]))
# plt.show()

mixer = Mixer(
    [tube_model.outlet_flow, tube_model_2.outlet_flow, tube_model_3.outlet_flow, tube_model_4.outlet_flow])
mixer.mass_balance()

for cmm in mixer.split_Out_flow.comp_dict:
    returnflow[cmm]=mixer.split_Out_flow.comp_dict[cmm]['mole_flow']

while np.sum(np.array([(tear_flow[c.name] - returnflow[c.name]) ** 2 for c in GlobalComponentManager.component_list]))>= 1e-12:
    tear_flow = returnflow
    for comps in flow_toR130.comp_dict:
        flow_toR130.comp_dict[comps]['mole_flow'] += tear_flow[comps]

    reactor = CstrSingleLiqPhase(198., 10132500., 2.8047, flow_toR130, reaction_set_1, properties_package,q_spec=6000)
    model1 = ReactorModel("R130", reactor, flow_toR130)

    reactor2 = CstrSingleLiqPhase(198., 10132500., 3, model1.outlet_flow, reaction_set_1, properties_package,q_spec=6000)
    model2 = ReactorModel("R130-2", reactor2, model1.outlet_flow)

    reactor3 = CstrSingleLiqPhase(198., 10132500., 3, model2.outlet_flow, reaction_set_1, properties_package,q_spec=6000)
    model3 = ReactorModel("R130-3", reactor3, model2.outlet_flow)

    reactor4 = CstrSingleLiqPhase(198., 10132500., 2.8047, model3.outlet_flow, reaction_set_1, properties_package,q_spec=6000)
    model4 = ReactorModel("R130-4", reactor4, model3.outlet_flow)

    solve_manager = SolverManager()
    solve_manager.add_model(model1)
    solve_manager.add_model(model2)
    solve_manager.add_model(model3)
    solve_manager.add_model(model4)
    # solve_manager.add_model(tube_model)
    result = solve_manager.solve_sequence()
    # PostProcess.process_results(result)
    # IB_flow = []
    # act_sp = []
    # for i in result:
    #     IB_flow.append(i['outflows']['IB'])
    #     act_sp.append(i['outflows']['p1_ion[1]'])
    # print(IB_flow)
    # plt.plot(act_sp)
    # plt.show()

    # print(pyo.value(model4.outlet_flow.comp_dict['IB']['mole_flow']))

    spliter_0 = Spliter(model4.outlet_flow, [0.001453043, 0.418959329, 0.274750957, 0.143179036, 0.161657635])
    spliter_0.mass_balance()
    # print(pyo.value(spliter_0.split_Out_flow_list[0].comp_dict['IP']['mole_flow']))
    tubeR = PFRSingleliqPhase(198, 10132500, 3.1, 0.18372636, spliter_0.split_Out_flow_list[1], reaction_set_1,
                              properties_package)
    tube_model = ReactorModel_PFR('tube_outer_1', tubeR, inlet_flow=spliter_0.split_Out_flow_list[1])
    solve_manager2 = SolverManager()
    solve_manager2.add_model(tube_model)
    res = solve_manager2.solve_sequence()
    PostProcess.process_results(res)

    tubeR_2 = PFRSingleliqPhase(198, 10132500, 3.1, 0.1836915, spliter_0.split_Out_flow_list[2], reaction_set_1,
                                properties_package)
    tube_model_2 = ReactorModel_PFR('tube_outer_2', tubeR_2, inlet_flow=spliter_0.split_Out_flow_list[2])
    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(tube_model_2.model, nfe=50, scheme='BACKWARD')
    solver2 = SolverFactory('ipopt')
    solver2.solve(tube_model_2.model)
    # print(tube_model_2.outlet_flow.comp_dict)
    tubeR_3 = PFRSingleliqPhase(198, 10132500, 3.1, 0.12245611, spliter_0.split_Out_flow_list[3], reaction_set_1,
                                properties_package)
    tube_model_3 = ReactorModel_PFR('tube_inner_1', tubeR_3, inlet_flow=spliter_0.split_Out_flow_list[3])
    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(tube_model_3.model, nfe=50, scheme='BACKWARD')
    solver2 = SolverFactory('ipopt')
    solver2.solve(tube_model_3.model)

    tubeR_4 = PFRSingleliqPhase(198, 10132500, 3.1, 0.12246347, spliter_0.split_Out_flow_list[4], reaction_set_1,
                                properties_package)
    tube_model_4 = ReactorModel_PFR('tube_inner_2', tubeR_4, inlet_flow=spliter_0.split_Out_flow_list[4])
    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(tube_model_4.model, nfe=50, scheme='BACKWARD')
    solver2 = SolverFactory('ipopt')
    solver2.solve(tube_model_4.model)

    # plt.plot(pyo.value(tube_model.model.F['IB', :]))
    # plt.plot(pyo.value(tube_model_2.model.F['IB', :]))
    # plt.plot(pyo.value(tube_model_3.model.F['IB', :]))
    # plt.plot(pyo.value(tube_model_4.model.F['IB', :]))
    # plt.show()

    mixer = Mixer(
        [tube_model.outlet_flow, tube_model_2.outlet_flow, tube_model_3.outlet_flow, tube_model_4.outlet_flow])
    mixer.mass_balance()

    for f in mixer.inlet_flow_list:
        print(f.comp_dict)

    for cmm in mixer.split_Out_flow.comp_dict:
        returnflow[cmm] = mixer.split_Out_flow.comp_dict[cmm]['mole_flow']

    print(np.sum(np.array([(tear_flow[c.name] - returnflow[c.name]) ** 2 for c in GlobalComponentManager.component_list])))
# comp_names = [c.name for c in GlobalComponentManager.component_list]
# init_values = {v: flow_toR130.comp_dict[v]['mole_flow'] + 1e-6 for v in flow_toR130.comp_dict}
# model.outflows = pyo.Var(comp_names, initialize=init_values, domain=pyo.NonNegativeReals)

print()
print(tear_flow)
print(returnflow)

# for c in flow_outR130.comp_dict:
#     flow_outR130.comp_dict[c]['mole_flow'] = model.outflows[c]
#
# model.eqs = pyo.ConstraintList()
#
# for eq in reactor.mass_balance(flow_outR130):
#     model.eqs.add(eq == 0.)
#
# solver = SolverFactory('ipopt')
# results = solver.solve(model, tee=True)
#
# if (results.solver.status == pyo.SolverStatus.ok) and (
#         results.solver.termination_condition == pyo.TerminationCondition.optimal):
#     print("Optimal solution found.")
#     for name in comp_names:
#         print(f"{name}: {pyo.value(model.outflows[name])}")
#     print("MWN: " + str(
#         (pyo.value(model.outflows['first_mom_live[1]'])+pyo.value(model.outflows['first_mom_dead[1]'])) / (pyo.value(model.outflows['zeroth_mom_dead[1]'])+pyo.value(model.outflows['zeroth_mom_dead[1]'])) * 56))
#     print("IB进料量：" + str(flow_toR130.comp_dict['IB']['mole_flow']))
#     print("转化率: " + str(
#         (flow_toR130.comp_dict['IB']['mole_flow'] - pyo.value(model.outflows['IB'])) / flow_toR130.comp_dict['IB'][
#             'mole_flow'] * 100))
# else:
#     print("Solver did not find an optimal solution.")
#
# print(str(pyo.value(model.outflows['first_mom_live[1]']) + pyo.value(model.outflows['first_mom_dead[1]'])))
