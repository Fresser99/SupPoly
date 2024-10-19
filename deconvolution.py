import pyomo.environ as pyo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyomo.opt import SolverFactory


class DeconvolutionHelper:

    def __init__(self, site_num=3, method='Gaussian'):
        self.site_num = site_num
        self.method = method
        self.MW = None
        self.Y = None

    def load_mwd(self, path):
        df = pd.read_csv(path)
        self.MW = np.array(df.iloc[:, 1].values)
        self.Y = np.array(df.iloc[:, -1].values)

    def plot_mwd(self):
        fig = plt.figure()
        plt.plot(np.log(self.MW), self.Y)
        plt.show()

    def deconvolution(self):

        model = pyo.ConcreteModel()
        param_name = []
        for cite in range(self.site_num):
            param_name.append('A' + str(cite))
            param_name.append('s' + str(cite))
            param_name.append('mu' + str(cite))
        initial_guess = {
            'A0': max(self.Y), 'A1': max(self.Y)/1.5, 'A2': max(self.Y)/2,
            'mu0': 0.,
            'mu1': 0.,
            'mu2': 0.,
            's0': 2, 's1': 2, 's2': 2
        }
        model.param = pyo.Var(param_name, initialize=initial_guess, domain=pyo.Reals)

        model.constraints = pyo.ConstraintList()
        for i in range(self.site_num):
            model.constraints.add(model.param['A0'] >= 0.0)
            model.constraints.add(model.param['A1'] >= 0.0)
            model.constraints.add(model.param['A2'] >= 0.0)
            model.constraints.add(model.param['A0'] <= 1.0)
            model.constraints.add(model.param['A1'] <= 1.0)
            model.constraints.add(model.param['A2'] <= 1.0)
            model.constraints.add(model.param['s' + str(i)] >= 0.2)
            model.constraints.add(model.param['mu' + str(i)] >= min(self.MW))
            model.constraints.add(model.param['mu' + str(i)] <= max(self.MW))

        def object_func(modelx):
            Y_cal = np.zeros([len(self.Y)])
            for c in range(self.site_num):
                Y_cal = Y_cal + self.gaussian(modelx, c)
            return np.sum((Y_cal - self.Y) ** 2)

        model.obj = pyo.Objective(rule=object_func, sense=pyo.minimize)
        solver = SolverFactory('ipopt')
        results = solver.solve(model, tee=True)
        peak1 = [pyo.value(model.param['A0']) * pyo.exp(
            -(self.MW[i] - pyo.value(model.param['mu0'])) ** 2 / (pyo.value(model.param['s0']) ** 2)) for i in
                 range(len(self.MW))]
        # peak1 = np.array(self.gaussian(model,0))
        # np.savetxt('fff.txt',peak1)
        peak2 = np.array([pyo.value(i) for i in self.gaussian(model, 1)])
        peak3 = np.array([pyo.value(i) for i in self.gaussian(model, 2)])
        # print([pyo.value(i) for i in self.gaussian(model, 0)])
        # plt.plot(peak1 + peak2 + peak3)
        # plt.plot(self.Y)
        plt.plot(peak1)
        plt.plot(peak2)
        plt.plot(peak3)
        # plt.plot(peak2)
        # plt.plot(peak3)
        plt.show()
        print([pyo.value(model.param['A' + str(i)]) for i in range(self.site_num)])
        print([pyo.value(model.param['mu' + str(i)]) for i in range(self.site_num)])
        print([pyo.value(model.param['s' + str(i)]) for i in range(self.site_num)])

    def gaussian(self, model, site_idx):

        a = [model.param['A' + str(site_idx)] * pyo.exp(
            -(self.MW[i] - model.param['mu' + str(site_idx)]) ** 2 / (
                    2 * model.param['s' + str(site_idx)] ** 2)) for i in range(len(self.MW))]
        return a


dc = DeconvolutionHelper()
dc.load_mwd('115504171830.csv')
dc.deconvolution()
