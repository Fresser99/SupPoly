import math

from pyomo import dae
from pyomo.common.enums import minimize
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
from solvermanage import *
import matplotlib.pyplot as plt
from OperationUnit import *
import pyomo.environ as pyo
from mwdrebuild import *

site_num = 3

IIR = Polymer([Segment('IB-R', 'IB-R-seg', CompType.segment), Segment('IP-R', 'IP-R-seg', CompType.segment)], 'IIR',
              'IIR', CompType.polymer)

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

# print(componentmanager.GlobalComponentManager.component_list)
for idx, i in enumerate(componentmanager.GlobalComponentManager.component_list):
    print(f'{idx}:' + i.name)

reaction_set_1 = ReactionSet()
kinetic_set = ['ka',
               'ki(1)_A', 'ki(2)', 'ki(3)',
               'kp(1)', 'kp(2)', 'kp(3)',
               'ktm(1)', 'ktm(2)', 'ktm(3)',
               'kd(1)', 'kd(2)', 'kd(3)']
param_estimation_problem = pyo.ConcreteModel()
param_estimation_problem.estimated_param = pyo.Var([k for k in kinetic_set],
                                                   domain=pyo.NonNegativeReals,
                                                   initialize={k: 0.000001 for k in kinetic_set})
kinetic = {'ka': param_estimation_problem.estimated_param['ka'],
           'ki(1)': param_estimation_problem.estimated_param['ki(1)'],
           'ki(2)': param_estimation_problem.estimated_param['ki(2)'],
           'ki(3)': param_estimation_problem.estimated_param['ki(3)'],
           'kp(1)': param_estimation_problem.estimated_param['kp(1)'],
           'kp(2)': param_estimation_problem.estimated_param['kp(2)'],
           'kp(3)': param_estimation_problem.estimated_param['kp(3)'],
           'ktm(1)': param_estimation_problem.estimated_param['ktm(1)'],
           'ktm(2)': param_estimation_problem.estimated_param['ktm(2)'],
           'ktm(3)': param_estimation_problem.estimated_param['ktm(3)'],
           'kd(1)': param_estimation_problem.estimated_param['kd(1)'],
           'kd(2)': param_estimation_problem.estimated_param['kd(2)'],
           'kd(3)': param_estimation_problem.estimated_param['kd(3)']}
# HCL + Ethy ——> p0 + ion-pair

# IB 守恒
# IB + ion-pair ——> P1[1]_A
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki_A(1)', 'value': kinetic['ki_A(1)']}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[2]_A
reaction_set_1.source_define(0, [0, 7], {'name': 'Ki_A(2)', 'value': kinetic['ki_A(2)']}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[3]_A
reaction_set_1.source_define(0, [0, 8], {'name': 'Ki_A(3)', 'value': kinetic['ki_A(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1]_A ——>  Y0[1]_A
reaction_set_1.source_define(0, [0, 15], {'name': 'Kp_AA(1)', 'value': kinetic['kp_AA(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2]_A ——>  Y0[2]_A
reaction_set_1.source_define(0, [0, 16], {'name': 'Kp_AA(2)', 'value': kinetic['kp_AA(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3]_A ——>  Y0[3]_A
reaction_set_1.source_define(0, [0, 17], {'name': 'Kp_AA(3)', 'value': kinetic['kp_AA(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1]_A ——> D0[1]_A+P1[1]_A
reaction_set_1.source_define(0, [0, 15], {'name': 'Ktm_A(1)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2]_A ——> D0[2]_A+P1[1]_A
reaction_set_1.source_define(0, [0, 16], {'name': 'Ktm_A(2)', 'value': kinetic['ktm_A(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3]_A ——> D0[3]_A+P1[3]_A
reaction_set_1.source_define(0, [0, 17], {'name': 'Ktm_A(3)', 'value': kinetic['ktm_A(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1]_B ——> Y0[1]_A
reaction_set_1.source_define(0, [0, 18], {'name': 'Kp_BA(1)', 'value': kinetic['kp_BA(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2]_B ——> Y0[1]_A
reaction_set_1.source_define(0, [0, 19], {'name': 'Kp_BA(2)', 'value': kinetic['kp_BA(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3]_B ——> Y0[1]_A
reaction_set_1.source_define(0, [0, 20], {'name': 'Kp_BA(3)', 'value': kinetic['kp_BA(3)']}, [1, 1], [1, 1], True)
# IB + Y0[1]_B ——> D0[1]_A+P1[1]_B
reaction_set_1.source_define(0, [0, 18], {'name': 'Ktm_BA(1)', 'value': kinetic['ktm_B(1)']}, [1, 1], [1, 1], True)
# IB + Y0[2]_B ——> D0[2]_A+P1[1]_B
reaction_set_1.source_define(0, [0, 19], {'name': 'Ktm_BA(2)', 'value': kinetic['ktm_B(2)']}, [1, 1], [1, 1], True)
# IB + Y0[3]_B ——> D0[3]_A+P1[3]_B
reaction_set_1.source_define(0, [0, 20], {'name': 'Ktm_BA(3)', 'value': kinetic['ktm_B(3)']}, [1, 1], [1, 1], True)

# IP
# Y0[1]_A + IP ——> Y0[1]_B
reaction_set_1.source_define(1, [1, 15], {'name': 'Kp_AB(1)', 'value': kinetic['kp_AB(1)']}, [1, 1], [1, 1], True)
# Y0[2]_A + IP ——> Y0[2]_B
reaction_set_1.source_define(1, [1, 16], {'name': 'Kp_AB(2)', 'value': kinetic['kp_AB(2)']}, [1, 1], [1, 1], True)
# Y0[3]_A + IP ——> Y0[3]_B
reaction_set_1.source_define(1, [1, 17], {'name': 'Kp_AB(3)', 'value': kinetic['kp_AB(3)']}, [1, 1], [1, 1], True)


# HCL
# HCL + EADC ——> ion-pair[1]
reaction_set_1.source_define(2, [2, 3], {'name': 'Ka(1)', 'value': kinetic['ka(1)']}, [1, 1], [1, 1], True)
# HCL + EADC ——> ion-pair[2]
reaction_set_1.source_define(2, [2, 3], {'name': 'Ka(2)', 'value': kinetic['ka(2)']}, [1, 1], [1, 1], True)
# HCL + EADC ——> ion-pair[3]
reaction_set_1.source_define(2, [2, 3], {'name': 'Ka(3)', 'value': kinetic['ka(3)']}, [1, 1], [1, 1], True)

# EADC
# EADC + HCL ——> ion-pair[1]
reaction_set_1.source_define(3, [2, 3], {'name': 'Ka(1)', 'value': kinetic['ka(1)']}, [1, 1], [1, 1], True)
# EADC + HCL ——> ion-pair[2]
reaction_set_1.source_define(3, [2, 3], {'name': 'Ka(2)', 'value': kinetic['ka(2)']}, [1, 1], [1, 1], True)
# EADC + HCL ——> ion-pair[3]
reaction_set_1.source_define(3, [2, 3], {'name': 'Ka(3)', 'value': kinetic['ka(3)']}, [1, 1], [1, 1], True)

# ion-pair
# HCL + EADC ——> ion-pair[1]
reaction_set_1.source_define(6, [2, 3], {'name': 'Ka(1)', 'value': kinetic['ka(1)']}, [1, 1], [1, 1], False)
# HCL + EADC ——> ion-pair[2]
reaction_set_1.source_define(7, [2, 3], {'name': 'Ka(2)', 'value': kinetic['ka(2)']}, [1, 1], [1, 1], False)
# HCL + EADC ——> ion-pair[3]
reaction_set_1.source_define(8, [2, 3], {'name': 'Ka(3)', 'value': kinetic['ka(3)']}, [1, 1], [1, 1], False)
# IB + ion-pair[1] ——> P1[1]_A
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki_A(1)', 'value': kinetic['ki_A(1)']}, [1, 1], [1, 1], True)
# IB + ion-pair[2] ——> P1[2]_A
reaction_set_1.source_define(7, [0, 7], {'name': 'Ki_A(2)', 'value': kinetic['ki_A(2)']}, [1, 1], [1, 1], True)
# IB + ion-pair[3] ——> P1[3]_A
reaction_set_1.source_define(8, [0, 8], {'name': 'Ki_A(3)', 'value': kinetic['ki_A(3)']}, [1, 1], [1, 1], True)

# p1_ion _A
# IB + ion-pair[1] ——> p1_ion[1]_A
reaction_set_1.source_define(9, [0, 6], {'name': 'Ki_A(1)', 'value': kinetic['ki_A(1)']}, [1, 1], [1, 1], False)
# IB + ion-pair[2] ——> p1_ion[2]_A
reaction_set_1.source_define(10, [0, 7], {'name': 'Ki_A(2)', 'value': kinetic['ki_A(2)']}, [1, 1], [1, 1], False)
# IB + ion-pair[3] ——> p1_ion[3]_A
reaction_set_1.source_define(11, [0, 8], {'name': 'Ki_A(3)', 'value': kinetic['ki_A(3)']}, [1, 1], [1, 1], False)
# IB + p1_ion[1]_A ——> p2[1]_A
reaction_set_1.source_define(9, [0, 9], {'name': 'Kp_AA(1)', 'value': kinetic['kp_AA(1)']}, [1, 1], [1, 1], True)
# IB + p1_ion[2]_A ——> p2[2]_A
reaction_set_1.source_define(10, [0, 10], {'name': 'Kp_AA(2)', 'value': kinetic['kp_AA(2)']}, [1, 1], [1, 1], True)
# IB + p1_ion[3]_A ——> p2[3]_A
reaction_set_1.source_define(11, [0, 11], {'name': 'Kp_AA(3)', 'value': kinetic['kp_AA(3)']}, [1, 1], [1, 1], True)
# IP + p1_ion[1]_A ——> p2[1]_B
reaction_set_1.source_define(9, [1, 9], {'name': 'Kp_AB(1)', 'value': kinetic['kp_AB(1)']}, [1, 1], [1, 1], True)
# IP + p1_ion[2]_A ——> p2[2]_B
reaction_set_1.source_define(10, [1, 10], {'name': 'Kp_AB(2)', 'value': kinetic['kp_AB(2)']}, [1, 1], [1, 1], True)
# IP + p1_ion[3]_A ——> p2[3]_B
reaction_set_1.source_define(11, [1, 11], {'name': 'Kp_AB(3)', 'value': kinetic['kp_AB(3)']}, [1, 1], [1, 1], True)
# Y0[1]_A + IB ——> Dn[1] + p1_ion[1]_A
reaction_set_1.source_define(9, [0, 15], {'name': 'Ktm_A(1)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], False)
# Y0[2]_A + IB ——> Dn[2] + p1_ion[2]_A
reaction_set_1.source_define(10, [0, 16], {'name': 'Ktm_A(2)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], False)
# Y0[3]_A + IB ——> Dn[2] + p1_ion[3]_A
reaction_set_1.source_define(11, [0, 17], {'name': 'Ktm_A(3)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], False)
# p1_ion[1]_A + IB ——> D1[1]: p1 monomer transfer
reaction_set_1.source_define(9, [0, 7], {'name': 'Ktm_A(1)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], True)
# p1_ion[2]_A + IB ——> D1[2]
reaction_set_1.source_define(10, [0, 8], {'name': 'Ktm_A(2)', 'value': kinetic['ktm_A(2)']}, [1, 1], [1, 1], True)
# p1_ion[3]_A + IB ——> D1[3]
reaction_set_1.source_define(11, [0, 9], {'name': 'Ktm_A(3)', 'value': kinetic['ktm_A(3)']}, [1, 1], [1, 1], True)
# p1_ion[1]_A + counter_ion ——> D1[1]
reaction_set_1.source_define(9, [9, 42], {'name': 'Kd_A(1)', 'value': kinetic['kd_A(1)']}, [1, 1], [1, 1], True)
# p1_ion[2]_A + counter_ion ——> D1[2]
reaction_set_1.source_define(10, [10, 42], {'name': 'Kd_A(2)', 'value': kinetic['kd_A(2)']}, [1, 1], [1, 1], True)
# p1_ion[3]_A + counter_ion ——> D1[3]
reaction_set_1.source_define(11, [11, 42], {'name': 'Kd_A(3)', 'value': kinetic['kd_A(3)']}, [1, 1], [1, 1], True)

# p1_ion_B
# IB + Y0[1]_B ——> p1_ion[1]_B
reaction_set_1.source_define(12, [0, 18], {'name': 'Ktm_B(1)', 'value': kinetic['ktm_B(1)']}, [1, 1], [1, 1], False)
# IB + Y0[2]_B ——> p1_ion[2]_B
reaction_set_1.source_define(13, [0, 19], {'name': 'Ktm_B(2)', 'value': kinetic['ktm_B(2)']}, [1, 1], [1, 1], False)
# IB + Y0[3]_B ——> p1_ion[3]_B
reaction_set_1.source_define(14, [0, 20], {'name': 'Ktm_B(3)', 'value': kinetic['ktm_B(3)']}, [1, 1], [1, 1], False)
# p1_ion_B[1] + IB ——> p2[1]_A
reaction_set_1.source_define(12, [0, 12], {'name': 'Kp_BA(1)', 'value': kinetic['kp_BA(1)']}, [1, 1], [1, 1], False)
# p1_ion_B[2] + IB ——> p2[2]_A
reaction_set_1.source_define(13, [0, 13], {'name': 'Kp_BA(2)', 'value': kinetic['kp_BA(2)']}, [1, 1], [1, 1], False)
# p1_ion_B[3] + IB ——> p2[3]_A
reaction_set_1.source_define(14, [0, 14], {'name': 'Kp_BA(3)', 'value': kinetic['kp_BA(3)']}, [1, 1], [1, 1], False)
# p1_ion_B[1] + IB ——> Dn[1]
reaction_set_1.source_define(12, [0, 12], {'name': 'Ktm_B(1)', 'value': kinetic['ktm_B(1)']}, [1, 1], [1, 1], True)
# p1_ion_B[2] + IB ——> Dn[2]
reaction_set_1.source_define(13, [0, 13], {'name': 'Ktm_B(2)', 'value': kinetic['ktm_B(2)']}, [1, 1], [1, 1], True)
# p1_ion_B[3] + IB ——> Dn[3]
reaction_set_1.source_define(14, [0, 14], {'name': 'Ktm_B(3)', 'value': kinetic['ktm_B(3)']}, [1, 1], [1, 1], True)
# p1_ion_B[1] +c_ion ——> Dn[1]
reaction_set_1.source_define(12, [12, 42], {'name': 'Kd_B(1)', 'value': kinetic['kd_B(1)']}, [1, 1], [1, 1], True)
# p1_ion_B[2] +c_ion ——> Dn[2]
reaction_set_1.source_define(13, [13, 42], {'name': 'Kd_B(2)', 'value': kinetic['kd_B(2)']}, [1, 1], [1, 1], True)
# p1_ion_B[3] +c_ion ——> Dn[3]
reaction_set_1.source_define(14, [14, 42], {'name': 'Kd_B(3)', 'value': kinetic['kd_B(3)']}, [1, 1], [1, 1], True)

# Y0_A
# ion-pair[1] + IB ——> P1_ion[1]_A ( Y0[1]_A )
reaction_set_1.source_define(15, [0, 6], {'name': 'Ki_A(1)', 'value': kinetic['ki_A(1)']}, [1, 1], [1, 1], False)
# ion-pair[2] + IB ——> P1_ion[2]_A ( Y0[2]_A )
reaction_set_1.source_define(16, [0, 7], {'name': 'Ki_A(2)', 'value': kinetic['ki_A(2)']}, [1, 1], [1, 1], False)
# ion-pair[3] + IB ——> P1_ion[3]_A ( Y0[3]_A )
reaction_set_1.source_define(17, [0, 8], {'name': 'Ki_A(3)', 'value': kinetic['ki_A(3)']}, [1, 1], [1, 1], False)

# # Y0[1]_A + IB ——> D0[1] + p1_ion[1]_A
# reaction_set_1.source_define(15, [0, 15], {'name': 'Ktm_A(1)', 'value': kinetic['ktm_A(1)']}, [1, 1], [1, 1], False)
# # Y0[2]_A + IB ——> D0[2] + pi_ion[2]_A
# reaction_set_1.source_define(16, [0, 16], {'name': 'Ktm_A(2)', 'value': kinetic['ktm_A(2)']}, [1, 1], [1, 1], False)
# # Y0[3]_A + IB ——> D0[3] + p1_ion[3]_A
# reaction_set_1.source_define(17, [0, 17], {'name': 'Ktm_A(3)', 'value': kinetic['ktm_A(3)']}, [1, 1], [1, 1], False)

# Y0[1]_A + IP ——> Y0[1]_B
reaction_set_1.source_define(15, [1, 15], {'name': 'Kp_AB(1)', 'value': kinetic['kp_AB(1)']}, [1, 1], [1, 1], True)
# Y0[2]_A + IP ——> Y0[2]_B
reaction_set_1.source_define(16, [1, 16], {'name': 'Kp_AB(2)', 'value': kinetic['kp_AB(2)']}, [1, 1], [1, 1], True)
# Y0[3]_A + IP ——> Y0[3]_B
reaction_set_1.source_define(17, [1, 17], {'name': 'Kp_AB(3)', 'value': kinetic['kp_AB(3)']}, [1, 1], [1, 1], True)

# Y0[1]_B + IB ——> Y0[1]_A
reaction_set_1.source_define(15, [0, 15], {'name': 'Kp_BA(1)', 'value': kinetic['kp_BA(1)']}, [1, 1], [1, 1], False)
# Y0[2]_B + IB ——> Y0[2]_A
reaction_set_1.source_define(16, [0, 16], {'name': 'Kp_BA(2)', 'value': kinetic['kp_BA(2)']}, [1, 1], [1, 1], False)
# Y0[3]_B + IB ——> Y0[3]_A
reaction_set_1.source_define(17, [0, 17], {'name': 'Kp_BA(3)', 'value': kinetic['kp_BA(3)']}, [1, 1], [1, 1], False)

# Yo[1]_A + c-ion ——> Dn[1]
reaction_set_1.source_define(15, [15, 42], {'name': 'Kd_A(1)', 'value': kinetic['kd_A(1)']}, [1, 1], [1, 1], True)
# Yo[2]_A + c-ion ——> Dn[2]
reaction_set_1.source_define(16, [16, 42], {'name': 'Kd_A(2)', 'value': kinetic['kd_A(1)']}, [1, 1], [1, 1], True)
# Yo[3]_A + c-ion ——> Dn[3]
reaction_set_1.source_define(17, [17, 42], {'name': 'Kd_A(3)', 'value': kinetic['kd_A(1)']}, [1, 1], [1, 1], True)


# Y0_B
# Y0[1]_A + IP ——> Y0[1]_B
reaction_set_1.source_define(18, [1, 15], {'name': 'Kp_AB(1)', 'value': kinetic['kp_AB(1)']}, [1, 1], [1, 1], False)
# Y0[2]_A + IP ——> Y0[2]_B
reaction_set_1.source_define(19, [1, 16], {'name': 'Kp_AB(2)', 'value': kinetic['kp_AB(2)']}, [1, 1], [1, 1], False)
# Y0[3]_A + IP ——> Y0[3]_B
reaction_set_1.source_define(20, [1, 17], {'name': 'Kp_AB(3)', 'value': kinetic['kp_AB(3)']}, [1, 1], [1, 1], False)

# Y0[1]_B + IB ——> Y0[1]_A
reaction_set_1.source_define(18, [0, 18], {'name': 'Kp_BA(1)', 'value': kinetic['kp_BA(1)']}, [1, 1], [1, 1], True)
# Y0[2]_B + IP ——> Y0[2]_A
reaction_set_1.source_define(19, [0, 19], {'name': 'Kp_BA(2)', 'value': kinetic['kp_BA(2)']}, [1, 1], [1, 1], True)
# Y0[3]_B + IP ——> Y0[3]_A
reaction_set_1.source_define(20, [0, 20], {'name': 'Kp_BA(3)', 'value': kinetic['kp_BA(3)']}, [1, 1], [1, 1], True)

# Y0[1]_B + c-ion ——> Dn[1]
reaction_set_1.source_define(18, [18, 42], {'name': 'Kd_B(1)', 'value': kinetic['kd_B(1)']}, [1, 1], [1, 1], True)
# Y0[2]_B + c-ion ——> Dn[2]
reaction_set_1.source_define(19, [19, 42], {'name': 'Kd_B(1)', 'value': kinetic['kd_B(1)']}, [1, 1], [1, 1], True)
# Y0[3]_B + c-ion ——> Dn[3]
reaction_set_1.source_define(20, [20, 42], {'name': 'Kd_B(1)', 'value': kinetic['kd_B(1)']}, [1, 1], [1, 1], True)




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
reaction_set_1.source_define(14, [14, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], True)
# Y1[3] + counter-ion ——> Dn[3] + ion-pair[3]
reaction_set_1.source_define(15, [15, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], True)

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
reaction_set_1.source_define(17, [17, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], True)
# Y2[3] +counter-ion+ ——> D2[3] + ion-pair[3]
reaction_set_1.source_define(18, [18, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], True)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(19, [0, 10], {'name': 'Ktm(1)', 'value': kinetic['ktm(1)']}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2]
reaction_set_1.source_define(20, [0, 11], {'name': 'Ktm(2)', 'value': kinetic['ktm(2)']}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(21, [0, 12], {'name': 'Ktm(3)', 'value': kinetic['ktm(3)']}, [1, 1], [1, 1], False)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(19, [10, 28], {'name': 'Kd(1)', 'value': kinetic['kd(1)']}, [1, 1], [1, 1], False)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(20, [11, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], False)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(21, [12, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], False)

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
reaction_set_1.source_define(23, [14, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], False)
# Y1[3] + counter_ion ——> D1[3] +ion-pair[3]
reaction_set_1.source_define(24, [15, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], False)

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
reaction_set_1.source_define(26, [17, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], False)
# Y2[3] + counter_ion ——> D2[3] +ion-pair[3]
reaction_set_1.source_define(27, [18, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], False)

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
reaction_set_1.source_define(28, [11, 28], {'name': 'Kd(2)', 'value': kinetic['kd(2)']}, [1, 1], [1, 1], True)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(28, [12, 28], {'name': 'Kd(3)', 'value': kinetic['kd(3)']}, [1, 1], [1, 1], True)
# reaction_set_1.preview_reaction_equations()

flow_toR130 = Flow(100, 103, 'R130_InLet')
flow_toR130.set_MassFlow_conventional(
    {'IB': 4370., 'IP': 0, 'HCL': 0.238, 'EADC': 3.8, 'HEXANE': 6.3, 'CH3CL': 8600})
component_index = [component for component in flow_toR130.comp_dict]
init_values = {component: 100. for
               component in flow_toR130.comp_dict}
init_values_r = {component: 1000. for
                 component in flow_toR130.comp_dict}
param_estimation_problem.global_model = pyo.Block()
# # global_model.tear_flow = pyo.Var(component_index,
#                                  initialize=init_values,
#                                  domain=pyo.NonNegativeReals)
param_estimation_problem.global_model.tear_flow = pyo.Var(component_index, domain=pyo.NonNegativeReals,
                                                          initialize=init_values_r)
param_estimation_problem.global_model.reactor_1 = pyo.Block()
param_estimation_problem.global_model.reactor_1.outflow = pyo.Var(component_index, initialize=init_values,
                                                                  domain=pyo.NonNegativeReals)
flow_to_reactor1 = Flow(198, 101325, 'to_reactor1')
for c in flow_to_reactor1.comp_dict:
    flow_to_reactor1.comp_dict[c]['mole_flow'] = flow_toR130.comp_dict[c]['mole_flow']

reactor1_model = CstrSingleLiqPhase(198, 101325, 10.26, flow_to_reactor1, reaction_set_1, properties_package)
Out_flow_reactor1 = Flow(reactor1_model.Temperature, reactor1_model.Pressure, 'OutFlow_reactor1')
for c in Out_flow_reactor1.comp_dict:
    Out_flow_reactor1.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.reactor_1.outflow[c]
param_estimation_problem.global_model.reactor_1.mass_balance = pyo.ConstraintList()
for eq in reactor1_model.mass_balance(Out_flow_reactor1):
    param_estimation_problem.global_model.reactor_1.mass_balance.add(eq == 0.)

# param_estimation_problem.global_model.reactor_2 = pyo.Block()
# param_estimation_problem.global_model.reactor_2.outflow = pyo.Var(component_index, initialize=init_values,
#                                                                   domain=pyo.NonNegativeReals)
# reactor2_model = CstrSingleLiqPhase(198, 101325, 3.3, Out_flow_reactor1, reaction_set_1, properties_package)
# Out_flow_reactor2 = Flow(reactor2_model.Temperature, reactor2_model.Pressure, 'OutFlow_reactor2')
# for c in Out_flow_reactor2.comp_dict:
#     Out_flow_reactor2.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.reactor_2.outflow[c]
# param_estimation_problem.global_model.reactor_2.mass_balance = pyo.ConstraintList()
# for eq in reactor2_model.mass_balance(Out_flow_reactor2):
#     param_estimation_problem.global_model.reactor_2.mass_balance.add(eq == 0.)
#
# param_estimation_problem.global_model.reactor_3 = pyo.Block()
# param_estimation_problem.global_model.reactor_3.outflow = pyo.Var(component_index, initialize=init_values,
#                                                                   domain=pyo.NonNegativeReals)
# reactor3_model = CstrSingleLiqPhase(198, 151325, 2, Out_flow_reactor2, reaction_set_1, properties_package)
# Out_flow_reactor3 = Flow(reactor3_model.Temperature, reactor3_model.Pressure, 'OutFlow_reactor3')
# for c in Out_flow_reactor3.comp_dict:
#     Out_flow_reactor3.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.reactor_3.outflow[c]
# param_estimation_problem.global_model.reactor_3.mass_balance = pyo.ConstraintList()
# for eq in reactor3_model.mass_balance(Out_flow_reactor3):
#     param_estimation_problem.global_model.reactor_3.mass_balance.add(eq == 0.)
# #
# param_estimation_problem.global_model.reactor_4 = pyo.Block()
# param_estimation_problem.global_model.reactor_4.outflow = pyo.Var(component_index, initialize=init_values,
#                                                                   domain=pyo.NonNegativeReals)
# reactor4_model = CstrSingleLiqPhase(198, 101325, 2.8, Out_flow_reactor3, reaction_set_1, properties_package)
# Out_flow_reactor4 = Flow(reactor4_model.Temperature, reactor4_model.Pressure, 'OutFlow_reactor4')
# for c in Out_flow_reactor4.comp_dict:
#     Out_flow_reactor4.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.reactor_4.outflow[c]
# param_estimation_problem.global_model.reactor_4.mass_balance = pyo.ConstraintList()
# for eq in reactor4_model.mass_balance(Out_flow_reactor4):
#     param_estimation_problem.global_model.reactor_4.mass_balance.add(eq == 0.)
#
# spliter = Spliter(Out_flow_reactor2, [0.001453043, 0.418959329, 0.274750957, 0.143179036, 0.161657635])
spliter = Spliter(Out_flow_reactor1, [1])
spliter.mass_balance()

# param_estimation_problem.global_model.pfr1 = pyo.Block()
#
# pfr1_model = PFRSingleliqPhase(t=198, p=101325, l=6.014, D=0.6, inflow=spliter.split_Out_flow_list[1],
#                                rx_set=reaction_set_1, prop=properties_package)
#
#
# def initialize_F(model, comp, z):
#     if z == 0:
#         return pfr1_model.Inflow.comp_dict[comp]['mole_flow']
#     else:
#         return pfr1_model.Inflow.comp_dict[comp]['mole_flow']
#
#
# param_estimation_problem.global_model.pfr1.z = dae.ContinuousSet(bounds=(0, pfr1_model.Length))
# param_estimation_problem.global_model.pfr1.F = pyo.Var(component_index, param_estimation_problem.global_model.pfr1.z,
#                                                        domain=pyo.NonNegativeReals,
#                                                        initialize=initialize_F)
# param_estimation_problem.global_model.pfr1.dFdz = dae.DerivativeVar(param_estimation_problem.global_model.pfr1.F,
#                                                                     wrt=param_estimation_problem.global_model.pfr1.z)
# Out_flow_pfr1 = Flow(101, 101, 'OutFlow_pfr1')
# param_estimation_problem.global_model.pfr1.initial_flow_constraint = pyo.ConstraintList()
# for c in Out_flow_pfr1.comp_dict:
#     Out_flow_pfr1.comp_dict[c] = param_estimation_problem.global_model.pfr1.F[
#         c, param_estimation_problem.global_model.pfr1.z.last()]
#     param_estimation_problem.global_model.pfr1.initial_flow_constraint.add(
#         spliter.split_Out_flow_list[1].comp_dict[c]['mole_flow'] == param_estimation_problem.global_model.pfr1.F[c, 0])
#
# param_estimation_problem.global_model.pfr1.dpn = pyo.Expression(param_estimation_problem.global_model.pfr1.z,
#                                                                 rule=pfr1_model.compute_dpn)
# param_estimation_problem.global_model.pfr1.V_flow = pyo.Expression(param_estimation_problem.global_model.pfr1.z,
#                                                                    rule=pfr1_model.volume_flow_rate_rule)
# param_estimation_problem.global_model.pfr1.C = pyo.Expression(component_index,
#                                                               param_estimation_problem.global_model.pfr1.z,
#                                                               rule=pfr1_model.concentration_rule)
# param_estimation_problem.global_model.pfr1.mass_balance = pyo.Constraint(component_index,
#                                                                          param_estimation_problem.global_model.pfr1.z,
#                                                                          rule=pfr1_model.mass_balance)
#
# discretizer = pyo.TransformationFactory('dae.finite_difference')
# discretizer.apply_to(param_estimation_problem.global_model.pfr1, nfe=20, scheme='BACKWARD')
#
# # param_estimation_problem.global_model.pfr2 = pyo.Block()
# #
# # pfr2_model = PFRSingleliqPhase(t=198, p=101325, l=6.014, D=0.1837, inflow=spliter.split_Out_flow_list[2],
# #                                rx_set=reaction_set_1, prop=properties_package)
# #
# #
# # def initialize_F2(model, comp, z):
# #     if z == 0:
# #         return pfr2_model.Inflow.comp_dict[comp]['mole_flow']
# #     else:
# #         return pfr2_model.Inflow.comp_dict[comp]['mole_flow']
# #
# #
# # param_estimation_problem.global_model.pfr2.z = dae.ContinuousSet(bounds=(0, pfr2_model.Length))
# # param_estimation_problem.global_model.pfr2.F = pyo.Var(component_index, param_estimation_problem.global_model.pfr2.z,
# #                                                        domain=pyo.NonNegativeReals,
# #                                                        initialize=initialize_F2)
# # param_estimation_problem.global_model.pfr2.dFdz = dae.DerivativeVar(param_estimation_problem.global_model.pfr2.F,
# #                                                                     wrt=param_estimation_problem.global_model.pfr2.z)
# # Out_flow_pfr2 = Flow(101, 101, 'OutFlow_pfr2')
# # param_estimation_problem.global_model.pfr2.initial_flow_constraint = pyo.ConstraintList()
# # for c in Out_flow_pfr2.comp_dict:
# #     Out_flow_pfr2.comp_dict[c] = param_estimation_problem.global_model.pfr2.F[
# #         c, param_estimation_problem.global_model.pfr2.z.last()]
# #     param_estimation_problem.global_model.pfr2.initial_flow_constraint.add(
# #         spliter.split_Out_flow_list[2].comp_dict[c]['mole_flow'] == param_estimation_problem.global_model.pfr2.F[c, 0])
# #
# # param_estimation_problem.global_model.pfr2.dpn = pyo.Expression(param_estimation_problem.global_model.pfr2.z,
# #                                                                 rule=pfr2_model.compute_dpn)
# # param_estimation_problem.global_model.pfr2.V_flow = pyo.Expression(param_estimation_problem.global_model.pfr2.z,
# #                                                                    rule=pfr2_model.volume_flow_rate_rule)
# # param_estimation_problem.global_model.pfr2.C = pyo.Expression(component_index,
# #                                                               param_estimation_problem.global_model.pfr2.z,
# #                                                               rule=pfr2_model.concentration_rule)
# # param_estimation_problem.global_model.pfr2.mass_balance = pyo.Constraint(component_index,
# #                                                                          param_estimation_problem.global_model.pfr2.z,
# #                                                                          rule=pfr2_model.mass_balance)
# #
# # discretizer2 = pyo.TransformationFactory('dae.finite_difference')
# # discretizer2.apply_to(param_estimation_problem.global_model.pfr2, nfe=20, scheme='BACKWARD')
# #
# # param_estimation_problem.global_model.pfr3 = pyo.Block()
# #
# # pfr3_model = PFRSingleliqPhase(t=198, p=101325, l=6.014, D= 0.12245611, inflow=spliter.split_Out_flow_list[3],
# #                                rx_set=reaction_set_1, prop=properties_package)
# #
# #
# # def initialize_F3(model, comp, z):
# #     if z == 0:
# #         return pfr3_model.Inflow.comp_dict[comp]['mole_flow']
# #     else:
# #         return pfr3_model.Inflow.comp_dict[comp]['mole_flow']
# #
# #
# # param_estimation_problem.global_model.pfr3.z = dae.ContinuousSet(bounds=(0, pfr3_model.Length))
# # param_estimation_problem.global_model.pfr3.F = pyo.Var(component_index, param_estimation_problem.global_model.pfr3.z,
# #                                                        domain=pyo.NonNegativeReals,
# #                                                        initialize=initialize_F3)
# # param_estimation_problem.global_model.pfr3.dFdz = dae.DerivativeVar(param_estimation_problem.global_model.pfr3.F,
# #                                                                     wrt=param_estimation_problem.global_model.pfr3.z)
# # Out_flow_pfr3 = Flow(101, 101, 'OutFlow_pfr3')
# # param_estimation_problem.global_model.pfr3.initial_flow_constraint = pyo.ConstraintList()
# # for c in Out_flow_pfr3.comp_dict:
# #     Out_flow_pfr3.comp_dict[c] = param_estimation_problem.global_model.pfr3.F[
# #         c, param_estimation_problem.global_model.pfr3.z.last()]
# #     param_estimation_problem.global_model.pfr3.initial_flow_constraint.add(
# #         spliter.split_Out_flow_list[3].comp_dict[c]['mole_flow'] == param_estimation_problem.global_model.pfr3.F[c, 0])
# #
# # param_estimation_problem.global_model.pfr3.dpn = pyo.Expression(param_estimation_problem.global_model.pfr3.z,
# #                                                                 rule=pfr3_model.compute_dpn)
# # param_estimation_problem.global_model.pfr3.V_flow = pyo.Expression(param_estimation_problem.global_model.pfr3.z,
# #                                                                    rule=pfr3_model.volume_flow_rate_rule)
# # param_estimation_problem.global_model.pfr3.C = pyo.Expression(component_index,
# #                                                               param_estimation_problem.global_model.pfr3.z,
# #                                                               rule=pfr3_model.concentration_rule)
# # param_estimation_problem.global_model.pfr3.mass_balance = pyo.Constraint(component_index,
# #                                                                          param_estimation_problem.global_model.pfr3.z,
# #                                                                          rule=pfr3_model.mass_balance)
# #
# # discretizer3 = pyo.TransformationFactory('dae.finite_difference')
# # discretizer3.apply_to(param_estimation_problem.global_model.pfr3, nfe=40, scheme='BACKWARD')
# #
# # param_estimation_problem.global_model.pfr4 = pyo.Block()
# #
# # pfr4_model = PFRSingleliqPhase(t=198, p=101325, l=6.014, D=0.12246347, inflow=spliter.split_Out_flow_list[4],
# #                                rx_set=reaction_set_1, prop=properties_package)
# #
# #
# # def initialize_F4(model, comp, z):
# #     if z == 0:
# #         return pfr4_model.Inflow.comp_dict[comp]['mole_flow']
# #     else:
# #         return pfr4_model.Inflow.comp_dict[comp]['mole_flow']
# #
# #
# # param_estimation_problem.global_model.pfr4.z = dae.ContinuousSet(bounds=(0, pfr4_model.Length))
# # param_estimation_problem.global_model.pfr4.F = pyo.Var(component_index, param_estimation_problem.global_model.pfr4.z,
# #                                                        domain=pyo.NonNegativeReals,
# #                                                        initialize=initialize_F4)
# # param_estimation_problem.global_model.pfr4.dFdz = dae.DerivativeVar(param_estimation_problem.global_model.pfr4.F,
# #                                                                     wrt=param_estimation_problem.global_model.pfr4.z)
# # Out_flow_pfr4 = Flow(101, 101, 'OutFlow_pfr4')
# # param_estimation_problem.global_model.pfr4.initial_flow_constraint = pyo.ConstraintList()
# # for c in Out_flow_pfr4.comp_dict:
# #     Out_flow_pfr4.comp_dict[c] = param_estimation_problem.global_model.pfr4.F[
# #         c, param_estimation_problem.global_model.pfr4.z.last()]
# #     param_estimation_problem.global_model.pfr4.initial_flow_constraint.add(
# #         spliter.split_Out_flow_list[4].comp_dict[c]['mole_flow'] == param_estimation_problem.global_model.pfr4.F[c, 0])
# #
# # param_estimation_problem.global_model.pfr4.dpn = pyo.Expression(param_estimation_problem.global_model.pfr4.z,
# #                                                                 rule=pfr2_model.compute_dpn)
# # param_estimation_problem.global_model.pfr4.V_flow = pyo.Expression(param_estimation_problem.global_model.pfr4.z,
# #                                                                    rule=pfr2_model.volume_flow_rate_rule)
# # param_estimation_problem.global_model.pfr4.C = pyo.Expression(component_index,
# #                                                               param_estimation_problem.global_model.pfr4.z,
# #                                                               rule=pfr4_model.concentration_rule)
# # param_estimation_problem.global_model.pfr4.mass_balance = pyo.Constraint(component_index,
# #                                                                          param_estimation_problem.global_model.pfr4.z,
# #                                                                          rule=pfr4_model.mass_balance)
# #
# # discretizer4 = pyo.TransformationFactory('dae.finite_difference')
# # discretizer4.apply_to(param_estimation_problem.global_model.pfr4, nfe=40, scheme='BACKWARD')
#
# return_flow = Flow(101, 101, 'return')
# # for c in return_flow.comp_dict:
# #     return_flow.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.pfr1.F[
# #                                                 c, param_estimation_problem.global_model.pfr1.z.last()] + \
# #                                             param_estimation_problem.global_model.pfr2.F[
# #                                                 c, param_estimation_problem.global_model.pfr2.z.last()] + \
# #                                             param_estimation_problem.global_model.pfr3.F[
# #                                                 c, param_estimation_problem.global_model.pfr3.z.last()] + \
# #                                             param_estimation_problem.global_model.pfr4.F[
# #                                                 c, param_estimation_problem.global_model.pfr4.z.last()]
#
# for c in return_flow.comp_dict:
#     return_flow.comp_dict[c]['mole_flow'] = param_estimation_problem.global_model.pfr1.F[
#         c, param_estimation_problem.global_model.pfr1.z.last()]
# param_estimation_problem.global_model.tear_constrain = pyo.ConstraintList()
#
# for c in component_index:
#     param_estimation_problem.global_model.tear_constrain.add(
#         param_estimation_problem.global_model.tear_flow[c] == return_flow.comp_dict[c]["mole_flow"])
# param_estimation_problem.global_model.out_let_mwn_1 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'first_mom_live[1]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'first_mom_dead[1]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_live[1]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_dead[1]'][
#                                                                                                 "mole_flow"]) * 56))
# param_estimation_problem.global_model.out_let_mww_1 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'second_mom_live[1]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'second_mom_dead[1]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_live[1]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_dead[1]'][
#                                                                                                 "mole_flow"]) * 56))
# param_estimation_problem.global_model.out_let_mwn_2 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'first_mom_live[2]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'first_mom_dead[2]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_live[2]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_dead[2]'][
#                                                                                                 "mole_flow"]) * 56))
# param_estimation_problem.global_model.out_let_mww_2 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'second_mom_live[2]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'second_mom_dead[2]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_live[2]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_dead[2]'][
#                                                                                                 "mole_flow"]) * 56))
#
# param_estimation_problem.global_model.out_let_mwn_3 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'first_mom_live[3]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'first_mom_dead[3]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_live[3]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'zeroth_mom_dead[3]'][
#                                                                                                 "mole_flow"]) * 56))
# param_estimation_problem.global_model.out_let_mww_3 = pyo.Expression(expr=pyo.log10((spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#                                                                                          'second_mom_live[3]'][
#                                                                                          "mole_flow"] +
#                                                                                      spliter.split_Out_flow_list[
#                                                                                          0].comp_dict[
#
#                                                                                          'second_mom_dead[3]'][
#                                                                                          "mole_flow"]) / (
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_live[3]'][
#                                                                                                 "mole_flow"] +
#                                                                                             spliter.split_Out_flow_list[
#                                                                                                 0].comp_dict[
#                                                                                                 'first_mom_dead[3]'][
#                                                                                                 "mole_flow"]) * 56))

# param_estimation_problem.global_model.M0 = pyo.Expression(expr=spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[1]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[1]']['mole_flow']+spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[2]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[2]']['mole_flow']+spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[3]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[3]']['mole_flow'])
#
# param_estimation_problem.global_model.M1 = pyo.Expression(expr=spliter.split_Out_flow_list[0].comp_dict['first_mom_live[1]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[1]']['mole_flow']+spliter.split_Out_flow_list[0].comp_dict['first_mom_live[2]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[2]']['mole_flow']+spliter.split_Out_flow_list[0].comp_dict['first_mom_live[3]']['mole_flow'] +
#     spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[3]']['mole_flow'])
#
# param_estimation_problem.global_model.out_let_pdi_1 = pyo.Expression(
#     expr=param_estimation_problem.global_model.out_let_mww_1 / param_estimation_problem.global_model.out_let_mwn_1)
# param_estimation_problem.global_model.out_let_pdi_2 = pyo.Expression(
#     expr=param_estimation_problem.global_model.out_let_mww_2 / param_estimation_problem.global_model.out_let_mwn_2)
# param_estimation_problem.global_model.out_let_pdi_3 = pyo.Expression(
#     expr=param_estimation_problem.global_model.out_let_mww_3 / param_estimation_problem.global_model.out_let_mwn_3)
# param_estimation_problem.Obj = pyo.ObjectiveList()
# param_estimation_problem.Obj = pyo.Objective(expr=((spliter.split_Out_flow_list[0].comp_dict['IB']["mole_flow"] / (
#         spliter.split_Out_flow_list[0].comp_dict['IB']["mole_flow"] +
#         spliter.split_Out_flow_list[0].comp_dict['CH3CL']["mole_flow"])*100-3.2) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_pdi_1 - 1.1375) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_mwn_1 - math.log10(
#                                                            967069)) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_pdi_2 - 2.423691467120772) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_mwn_2 - math.log10(
#                                                            79585.50120483543)) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_pdi_3 - 1.524951272223372) ** 2 + (
#                                                                param_estimation_problem.global_model.out_let_mwn_3 - math.log10(
#                                                            335145.4219643992)) ** 2 +
#                                                    0.89*(pyo.log10(param_estimation_problem.global_model.M1/param_estimation_problem.global_model.M0*56)-math.log10(96000))**2), sense=minimize)

# param_estimation_problem.Obj.add(pyo.Objective(expr=((param_estimation_problem.global_model.reactor_4.outflow['first_mom_live[1]']+param_estimation_problem.global_model.reactor_4.outflow['first_mom_dead[1]'])/(param_estimation_problem.global_model.reactor_4.outflow['zeroth_mom_live[1]']+param_estimation_problem.global_model.reactor_4.outflow['zeroth_mom_dead[1]'])*56-96000)**2,sense=pyo.minimize))
solver = SolverFactory('ipopt')
solver.solve(param_estimation_problem, tee=True)

for c in param_estimation_problem.estimated_param.keys():
    print(f"{c}:" + str(pyo.value(param_estimation_problem.estimated_param[c])))
for c in component_index:
    print(f"{c}：" + str(pyo.value(spliter.split_Out_flow_list[0].comp_dict[c]['mole_flow'])))
# print(str(pyo.value(param_estimation_problem.global_model.reactor_1.outflow['IB'])))
# print(str(pyo.value(param_estimation_problem.global_model.reactor_2.outflow['IB'])))
# print(str(pyo.value(param_estimation_problem.global_model.reactor_3.outflow['IB'])))
print("反应器4出口IB：" + str(pyo.value(spliter.split_Out_flow_list[0].comp_dict['IB']["mole_flow"])))
print("反应器4出口CH3CL：" + str(pyo.value(spliter.split_Out_flow_list[0].comp_dict['CH3CL']["mole_flow"])))
print("反应器4出口mwn：" + str((pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['first_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[1]']['mole_flow'])) / (pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[1]']['mole_flow'])) * 56))
print("反应器4出口mww：" + str((pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['second_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['second_mom_dead[1]']['mole_flow'])) / (pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['first_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[1]']['mole_flow'])) * 56))

M0 = (pyo.value(spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[1]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[2]']['mole_flow']) + pyo.value(
             spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[2]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_live[3]']['mole_flow']) + pyo.value(
             spliter.split_Out_flow_list[0].comp_dict['zeroth_mom_dead[3]']['mole_flow']))

M1 = (pyo.value(spliter.split_Out_flow_list[0].comp_dict['first_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[1]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['first_mom_live[2]']['mole_flow']) + pyo.value(
             spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[2]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['first_mom_live[3]']['mole_flow']) + pyo.value(
             spliter.split_Out_flow_list[0].comp_dict['first_mom_dead[3]']['mole_flow']))

M2 = (pyo.value(spliter.split_Out_flow_list[0].comp_dict['second_mom_live[1]']['mole_flow']) + pyo.value(
    spliter.split_Out_flow_list[0].comp_dict['second_mom_dead[1]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['second_mom_live[2]']['mole_flow']) + pyo.value(
             spliter.split_Out_flow_list[0].comp_dict['second_mom_dead[2]']['mole_flow'])) + (
                 pyo.value(spliter.split_Out_flow_list[0].comp_dict['second_mom_live[3]']['mole_flow']) + pyo.value(

             spliter.split_Out_flow_list[0].comp_dict['second_mom_dead[3]']['mole_flow']))
mwd = PolymerMWD(M0, M1, M2)

# 打印分布参数
params = mwd.get_distribution_parameters()
print("\n分布参数:")
for key, value in params.items():
    print(f"{key}: {value:.2e}" if isinstance(value, float) and value > 1000 else f"{key}: {value:.4f}")
mwd.plot_distribution()
