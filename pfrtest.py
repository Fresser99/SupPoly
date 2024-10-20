import pyomo.dae as dae
import pyomo.environ as pyo
from componentmanager import *
model = pyo.ConcreteModel()
model.z = dae.ContinuousSet(bounds=(0, 4))
comp_names = [c.name for c in GlobalComponentManager.component_list]
model.comps = pyo.Set(initialize=comp_names)
model.F = pyo.Var(model.comps, model.z, domain=pyo.NonNegativeReals)
model.dFdz = dae.DerivativeVar(model.F, wrt=model.z)
print(model.F['IB',model.z.first()])
