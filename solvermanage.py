from typing import List, Dict, Any

from pyomo import dae
from pyomo.opt import SolverFactory

from flow import *
import pyomo.environ as pyo


class ReactorModel:

    def __init__(self, name, model, inlet_flow):
        self.name = name
        self.inlet_flow = inlet_flow
        self.outlet_flow = Flow(t=0, p=101, name=f"{name}_out")
        self.reactor = model
        self.model = pyo.ConcreteModel()
        self.setup_model()

    def setup_model(self):
        comp_names = [c.name for c in GlobalComponentManager.component_list]
        init_values = {v: self.inlet_flow.comp_dict[v]['mole_flow'] + 1e-6 for v in self.inlet_flow.comp_dict}
        self.model.outflows = pyo.Var(comp_names, initialize=init_values, domain=pyo.NonNegativeReals)

        for c in self.outlet_flow.comp_dict:
            self.outlet_flow.comp_dict[c]['mole_flow'] = self.model.outflows[c]

        self.model.eqs = pyo.ConstraintList()
        for eq in self.reactor.mass_balance(self.outlet_flow):
            self.model.eqs.add(eq == 0.)

    def solve(self, solver):
        results = solver.solve(self.model, tee=True)
        if (results.solver.status == pyo.SolverStatus.ok) and (
                results.solver.termination_condition == pyo.TerminationCondition.optimal):
            return {name: pyo.value(self.model.outflows[name]) for name in self.model.outflows}
        else:
            raise RuntimeError(f"Solver failed to find an optimal solution for {self.name}")


class ReactorModel_PFR:

    def __init__(self, name, model, inlet_flow):
        self.name = name
        self.inlet_flow = inlet_flow
        self.outlet_flow = Flow(t=0, p=101, name=f"{name}_out")
        self.reactor = model
        self.model = pyo.ConcreteModel()
        self.setup_model()

    def setup_model(self):
        self.model.z = dae.ContinuousSet(bounds=(0, self.reactor.Length))
        comp_names = [c.name for c in GlobalComponentManager.component_list]
        self.model.comps = pyo.Set(initialize=comp_names)

        def initialize_F(model, comp, z):
            if z == 0:
                return self.inlet_flow.comp_dict[comp]['mole_flow']+1e-32
            else:
                return self.inlet_flow.comp_dict[comp]['mole_flow']+1e-32
        self.model.F = pyo.Var(self.model.comps, self.model.z, domain=pyo.NonNegativeReals,initialize=initialize_F)
        self.model.dFdz = dae.DerivativeVar(self.model.F, wrt=self.model.z)

        self.model.initial_flow_constraint = pyo.ConstraintList()

        for comp in self.model.comps:

            self.model.initial_flow_constraint.add(
                self.model.F[comp, 0] == self.inlet_flow.comp_dict[comp]['mole_flow'])

        def compute_dpn(model, z):
            polymer_zeroth_mole_flow = sum(
                self.model.F[c, z] if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] == 0 else 0.0 for c in
                self.inlet_flow.comp_dict)

            polymer_first_mole_flow = sum(
                self.model.F[c, z] if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] == 1 else 0.0 for c in
                self.inlet_flow.comp_dict)

            mole_flow_zeroth = sum(
                self.model.F[c, z] if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] in [0, -2] else 0.0 for c in
                self.inlet_flow.comp_dict)

            mole_flow_first = sum(
                self.model.F[c, z] if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] in [1, -2] else 0.0 for c in
                self.inlet_flow.comp_dict)

            return polymer_first_mole_flow / polymer_zeroth_mole_flow

        self.model.dpn = pyo.Expression(self.model.z, rule=compute_dpn)

        def volume_flow_rate_rule(model, z):

            polymer_zeroth_mole_flow = np.array(
                [pyo.value(self.model.F[c, z]) if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] == 0 else 0.0
                 for c in
                 self.inlet_flow.comp_dict])

            mole_flow_zeroth = np.array([pyo.value(self.model.F[c, z]) if self.inlet_flow.comp_dict[c][
                                                                              'polymer_flow_momentum'] in [0,
                                                                                                           -2] else 0.0
                                         for c in
                                         self.inlet_flow.comp_dict])

            mole_flow_first = [
                pyo.value(self.model.F[c, z]) if self.inlet_flow.comp_dict[c]['polymer_flow_momentum'] in [1,
                                                                                                           -2] else 0.0
                for c in
                self.inlet_flow.comp_dict]

            c_idx_list = []
            for c_idx, comp in enumerate(GlobalComponentManager.component_list):
                if type(comp) is Component:
                    c_idx_list.append(c_idx)

            mole_flow_properties = np.append(mole_flow_zeroth[c_idx_list], np.sum(polymer_zeroth_mole_flow))
            mole_frac_properties = mole_flow_properties / np.sum(mole_flow_properties)

            pc_ftr_polymer = self.reactor.PropertiesMethod.param.r[-1]

            self.reactor.PropertiesMethod.param.m[-1] = pc_ftr_polymer * pyo.value(self.model.dpn[z]) * \
                                                        self.reactor.PropertiesMethod.param.MW[-1]


            vm_liq = self.reactor.PropertiesMethod.calculate_molar_density_mixture(self.reactor.Temperature,
                                                                                   self.reactor.Pressure,
                                                                                   self.reactor.PropertiesMethod.param,
                                                                                   mole_frac_properties,
                                                                                   mole_frac_properties[-1],
                                                                                   self.model.dpn[z])/1000
            print(pyo.value(vm_liq))

            return vm_liq

        self.model.V_flow = pyo.Expression(self.model.z, rule=volume_flow_rate_rule)

        def concentration_rule(model, comp, z):
            return model.F[comp, z] / model.V_flow[z]

        self.model.C = pyo.Expression(self.model.comps, self.model.z, rule=concentration_rule)
        self.model.mass_balance = pyo.Constraint(self.model.comps, self.model.z, rule=self.reactor.mass_balance)

    def solve(self, solver):
        discretizer = pyo.TransformationFactory('dae.finite_difference')
        discretizer.apply_to(self.model, nfe=30, scheme='BACKWARD')
        results = solver.solve(self.model, tee=True)
        if (results.solver.status == pyo.SolverStatus.ok) and (
                results.solver.termination_condition == pyo.TerminationCondition.optimal):
            return {name: pyo.value(self.model.F[name, -1]) for name in self.model.comps}
        else:
            raise RuntimeError(f"Solver failed to find an optimal solution for {self.name}")


class SolverManager:
    def __init__(self):
        self.model_sequence: List[ReactorModel] = []
        self.solver = SolverFactory('ipopt')

    def add_model(self, model: ReactorModel):
        self.model_sequence.append(model)

    def solve_sequence(self) -> List[Dict[str, Any]]:
        results = []
        for i, model in enumerate(self.model_sequence):
            if i > 0:
                # Update inlet flow of current model with outlet flow of previous model
                for c in model.inlet_flow.comp_dict:
                    model.inlet_flow.comp_dict[c]['mole_flow'] = results[-1]['outflows'][c]
                model.setup_model()  # Reinitialize the model with new inlet conditions

            model_results = model.solve(self.solver)
            results.append({'name': model.name, 'outflows': model_results, 'inlet_flow': model.inlet_flow})
        return results


class PostProcess:
    @staticmethod
    def calculate_mwn(result: Dict[str, float]) -> float:
        return ((result['first_mom_live[1]'] + result['first_mom_dead[1]']) /
                (result['zeroth_mom_dead[1]'] + result['zeroth_mom_dead[1]']) * 56)

    @staticmethod
    def calculate_conversion(inlet_flow: Dict[str, float], outlet_flow: Dict[str, float],
                             component: str = 'IB') -> float:

        initial = inlet_flow[component]['mole_flow']
        final = outlet_flow[component]
        return (initial - final) / initial * 100

    @staticmethod
    def process_results(results: List[Dict[str, Any]]):
        for result in results:
            print(f"\n反应器 {result['name']} 的结果：")
            print("出料组成：")
            for name, value in result['outflows'].items():
                print(f"{name}: {value}")

            mwn = PostProcess.calculate_mwn(result['outflows'])
            print(f"数均分子量（MWN）: {mwn:.2f}")

            conversion = PostProcess.calculate_conversion(result['inlet_flow'].comp_dict, result['outflows'])
            print(f"IB转化率: {conversion:.2f}%")

            # 计算总体转化率
        initial_IB = results[0]['inlet_flow'].comp_dict['IB']['mole_flow']
        final_IB = results[-1]['outflows']['IB']
        overall_conversion = (initial_IB - final_IB) / initial_IB * 100
        print(f"\n总体IB转化率: {overall_conversion:.2f}%")

        # 最终的数均分子量（MWN）
        final_mwn = PostProcess.calculate_mwn(results[-1]['outflows'])
        print(f"最终数均分子量（MWN）: {final_mwn:.2f}")
