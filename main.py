import componentmanager
from component import *
from reactor import CstrSingleLiqPhase
from reactions import *
from flow import *
from proptiesmethod import *
import numpy as np

site_num = 3

IIR = Polymer([Segment('IB-R', 'IB-R-seg', CompType.segment)], 'IIR', 'IIR', CompType.polymer)

component_list = [Component("IB", '115-11-7', CompType.conventional), Component("IP", '78-79-5', CompType.conventional),
                  Component("HCL", '7647-01-0', CompType.conventional),
                  Component("EADC", '563-43-9', CompType.conventional),
                  Component("HEXANE", '110-54-3', CompType.conventional),
                  Component("CH3CL", '110-54-3', CompType.conventional), IIR]

componentmanager.GlobalComponentManager.component_list_gen(component_list, 'CATION', site_num)

print(len(componentmanager.GlobalComponentManager.component_list))
for c in GlobalComponentManager.component_list:
    if type(c) is Component:
        print(c.Formular)
    elif type(c) is SpecificComponent:
        print(c.name)
reaction_set_1 = ReactionSet()

# IB + ion-pair ——> P1[1] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[2] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + ion-pair ——> P1[3] + counter-ion
reaction_set_1.source_define(0, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[1] ——>  Y0[1]
reaction_set_1.source_define(0, [0, 10], {'name': 'Kp(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[2] ——>  Y0[2]
reaction_set_1.source_define(0, [0, 11], {'name': 'Kp(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[3] ——>  Y0[3]
reaction_set_1.source_define(0, [0, 12], {'name': 'Kp(3)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[1] ——> D0[1]+P1[1]
reaction_set_1.source_define(0, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[2] ——> D0[2]+P1[1]
reaction_set_1.source_define(0, [0, 11], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + Y0[3] ——> D0[3]+P1[3]
reaction_set_1.source_define(0, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# HCL + EADC ——> ion-pair
reaction_set_1.source_define(2, [2, 3], {'name': 'Ka', 'value': 0.0}, [1, 1], [1, 1], True)
# EADC + HCL ——> ion-pair
reaction_set_1.source_define(3, [2, 3], {'name': 'Ka', 'value': 0.0}, [1, 1], [1, 1], True)

# HCL + EADC ——> ion-pair
reaction_set_1.source_define(6, [2, 3], {'name': 'Ka', 'value': 0.0}, [1, 1], [1, 1], False)

# IB + ion-pair ——> P1[1] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], True)

# IB + ion-pair ——> P1[2] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], True)

# IB + ion-pair ——> P1[3] + counter-ion
reaction_set_1.source_define(6, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + conuter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(6, [10, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + conuter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(6, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + conuter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(6, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# IB + ion-pair ——> p1_ion[1] + counter-ion
reaction_set_1.source_define(7, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# IB + ion-pair ——> p1_ion[2] + counter-ion
reaction_set_1.source_define(8, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# IB + ion-pair ——> p1_ion[3] + counter-ion
reaction_set_1.source_define(9, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# IB + p1_ion[1] ——> p2[1]
reaction_set_1.source_define(7, [0, 7], {'name': 'Kp(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + p1_ion[2] ——> p2[2]
reaction_set_1.source_define(8, [0, 8], {'name': 'Kp(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# IB + p1_ion[3] ——> p2[3]
reaction_set_1.source_define(9, [0, 9], {'name': 'Kp(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# p1_ion[1] + IB ——> D1[1]: p1 monomer transfer
reaction_set_1.source_define(7, [0, 7], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> D1[2]
reaction_set_1.source_define(8, [0, 8], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> D1[3]
reaction_set_1.source_define(9, [0, 9], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# p1_ion[1] + counter_ion ——> D1[1] + ion-pair
reaction_set_1.source_define(7, [7, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[2] + counter_ion ——> D1[2] + ion-pair
reaction_set_1.source_define(8, [8, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + counter_ion ——> D1[3] + ion-pair
reaction_set_1.source_define(9, [9, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(7, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2]
reaction_set_1.source_define(8, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(9, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# ion-pair + IB ——> P1_ion[1] ( Y0[1] )
reaction_set_1.source_define(10, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> P1_ion[2] ( Y0[2] )
reaction_set_1.source_define(11, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> P1_ion[3] ( Y0[3] )
reaction_set_1.source_define(12, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(10, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + pi_ion[2]
reaction_set_1.source_define(11, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(12, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(10, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[2] + IB ——> D0[2] + pi_ion[2]
reaction_set_1.source_define(11, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[3] + IB ——> D0[3] + pi_ion[3]
reaction_set_1.source_define(12, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(10, [10, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(11, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(12, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] (1*Y0[1]=Y1[1] n=1)
reaction_set_1.source_define(13, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] (1*Y0[2]=Y1[2] n=1)
reaction_set_1.source_define(14, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] (1*Y0[3]=Y1[3] n=1)
reaction_set_1.source_define(15, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1](1*Y0[1]=Y1[1] n=1)
reaction_set_1.source_define(13, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2](1*Y0[2]=Y1[2] n=1)
reaction_set_1.source_define(14, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3](1*Y0[3]=Y1[3] n=1)
reaction_set_1.source_define(15, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y1[1] + IB ——> Y1[1]
reaction_set_1.source_define(13, [0, 10], {'name': 'Kp(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[2] + IB ——> Y1[2]
reaction_set_1.source_define(14, [0, 11], {'name': 'Kp(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[3] + IB ——> Y1[3]
reaction_set_1.source_define(15, [0, 12], {'name': 'Kp(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y1[1] + IB ——> Dn[1] + p1_ion[1]
reaction_set_1.source_define(13, [0, 13], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y1[2] + IB ——> Dn[2] + p1_ion[2]
reaction_set_1.source_define(14, [0, 14], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y1[3] + IB ——> Dn[3] + p1_ion[3]
reaction_set_1.source_define(15, [0, 15], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y1[1] + counter-ion ——> Dn[1] + ion-pair[1]
reaction_set_1.source_define(13, [13, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y1[2] + counter-ion ——> Dn[2] + ion-pair[2]
reaction_set_1.source_define(14, [14, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y1[3] + counter-ion ——> Dn[3] + ion-pair[3]
reaction_set_1.source_define(15, [15, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] (1^2*Y0[1]=Y2[1] n=1)
reaction_set_1.source_define(16, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] (1^2*Y0[2]=Y2[2] n=1)
reaction_set_1.source_define(17, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] (1^2*Y0[3]=Y2[3] n=1)
reaction_set_1.source_define(18, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + IB ——> D0[1] + p1_ion[1](1^2*Y0[1]=Y2[1] n=1)
reaction_set_1.source_define(16, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2](1^2*Y0[1]=Y2[2] n=1)
reaction_set_1.source_define(17, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3](1^2*Y0[1]=Y2[3] n=1)
reaction_set_1.source_define(18, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# 2Y1[1] + IB ——> Y2[1]
reaction_set_1.source_define(16, [0, 13], {'name': 'Kp(1)', 'value': 0.0}, [1, 2], [1, 1], False)
# 2Y1[2] + IB ——> Y2[2]
reaction_set_1.source_define(17, [0, 14], {'name': 'Kp(2)', 'value': 0.0}, [1, 2], [1, 1], False)
# 2Y1[3] + IB ——> Y2[3]
reaction_set_1.source_define(18, [0, 15], {'name': 'Kp(3)', 'value': 0.0}, [1, 2], [1, 1], False)

# Y0[1] + IB ——> Y2[1]
reaction_set_1.source_define(16, [0, 10], {'name': 'Kp(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> Y2[2]
reaction_set_1.source_define(17, [0, 11], {'name': 'Kp(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> Y2[3]
reaction_set_1.source_define(18, [0, 12], {'name': 'Kp(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y2[1] + IB ——> D2[1] + p1_ion[1]
reaction_set_1.source_define(16, [0, 16], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y2[2] + IB ——> D2[1] + p1_ion[2]
reaction_set_1.source_define(17, [0, 17], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y2[3] + IB ——> D2[3] + p1_ion[3]
reaction_set_1.source_define(18, [0, 18], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y2[1] +counter-ion+ ——> D2[1] + ion-pair[1]
reaction_set_1.source_define(16, [16, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y2[2] +counter-ion+ ——> D2[2] + ion-pair[2]
reaction_set_1.source_define(17, [17, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y2[3] +counter-ion+ ——> D2[3] + ion-pair[3]
reaction_set_1.source_define(18, [18, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y0[1] + IB ——> D0[1] + p1_ion[1]
reaction_set_1.source_define(19, [0, 10], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + IB ——> D0[2] + p1_ion[2]
reaction_set_1.source_define(20, [0, 11], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + IB ——> D0[3] + p1_ion[3]
reaction_set_1.source_define(21, [0, 12], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(19, [10, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(20, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(21, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(19, [0, 7], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(20, [0, 8], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(21, [0, 9], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(22, [0, 13], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(23, [0, 14], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[1] + IB ——> D1[1] +p1_ion[1]
reaction_set_1.source_define(24, [0, 15], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y1[1] + counter_ion ——> D1[1] +ion-pair[1]
reaction_set_1.source_define(22, [13, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[2] + counter_ion ——> D1[2] +ion-pair[2]
reaction_set_1.source_define(23, [14, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y1[3] + counter_ion ——> D1[3] +ion-pair[3]
reaction_set_1.source_define(24, [15, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(22, [0, 7], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(23, [0, 8], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(24, [0, 9], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# Y2[1] + IB ——> D2[1] + p1_ion[1]
reaction_set_1.source_define(25, [0, 22], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y2[2] + IB ——> D2[2] + p1_ion[2]
reaction_set_1.source_define(26, [0, 23], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y2[3] + IB ——> D2[3] + p1_ion[3]
reaction_set_1.source_define(27, [0, 24], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y2[1] + counter_ion ——> D2[1] +ion-pair[1]
reaction_set_1.source_define(25, [16, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y2[2] + counter_ion ——> D2[2] +ion-pair[2]
reaction_set_1.source_define(26, [17, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# Y2[3] + counter_ion ——> D2[3] +ion-pair[3]
reaction_set_1.source_define(27, [18, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# p1_ion[1] + IB ——> IB [1]
reaction_set_1.source_define(25, [0, 7], {'name': 'Ktm(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[2] + IB ——> IB [2]
reaction_set_1.source_define(26, [0, 8], {'name': 'Ktm(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# p1_ion[3] + IB ——> IB [3]
reaction_set_1.source_define(27, [0, 9], {'name': 'Ktm(3)', 'value': 0.0}, [1, 1], [1, 1], True)

# ion-pair + IB ——> p1_ion[1] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(1)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[2] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(2)', 'value': 0.0}, [1, 1], [1, 1], False)
# ion-pair + IB ——> p1_ion[3] + counter-ion
reaction_set_1.source_define(28, [0, 6], {'name': 'Ki(3)', 'value': 0.0}, [1, 1], [1, 1], False)

# Y0[1] + counter-ion ——> D0[1] + ion-pair
reaction_set_1.source_define(28, [10, 28], {'name': 'Kd(1)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[2] + counter-ion ——> D0[2] + ion-pair
reaction_set_1.source_define(28, [11, 28], {'name': 'Kd(2)', 'value': 0.0}, [1, 1], [1, 1], True)
# Y0[3] + counter-ion ——> D0[3] + ion-pair
reaction_set_1.source_define(28, [12, 28], {'name': 'Kd(3)', 'value': 0.0}, [1, 1], [1, 1], True)

flow_toR130 = Flow(100, 103, 'to_R130')

properties_method = UserMethod()
reactor = CstrSingleLiqPhase(100., 100., 30, flow_toR130, reaction_set_1, properties_method)
print(reaction_set_1.source_dict)
reaction_set_1.preview_reaction_equations()
