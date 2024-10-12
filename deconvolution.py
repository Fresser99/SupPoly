import pyomo.environ as pyo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class DeconvolutionHelper:

    def __init__(self, site_num=3, method='Gaussian'):
        self.site_num = site_num
        self.method = method
        self.MW = None
        self.Y = None

    def load_mwd(self, path):
        df = pd.read_csv(path)
        self.MW = np.array(df.iloc[:, 0].values)
        self.Y = np.array(df.iloc[:, -1].values)

    def plot_mwd(self):
        fig = plt.figure()
        plt.plot(np.log(self.MW), self.Y)
        plt.show()

    def deconvolution(self):
        pass




dc = DeconvolutionHelper()
dc.load_mwd('115504171830.csv')
dc.plot_mwd()
