from abc import ABC, abstractmethod
import numpy as np


class PropertiesMethod(ABC):

    @abstractmethod
    def calculate_molar_density_mixture(self, *args, **kwargs):
        pass

    @abstractmethod
    def calculate_molar_density_pure(self, *args, **kwargs):
        pass

    @abstractmethod
    def calculate_polymer_density(self, *args, **kwargs):
        pass


class UserMethod(PropertiesMethod):

    def calculate_molar_density_mixture(self, mole_frac, temperature, pressure, param_list):
        R = 8.314
        Tr = temperature / param_list["Tc"]
        m = 0.37464 + 1.54226 * param_list["Omega"] - 0.26992 * param_list["Omega"] ** 2
        alpha = (1 + m * (1 - np.sqrt(Tr))) ** 2

        a = 0.45724 * R ** 2 * param_list["Tc"] ** 2 / param_list["Pc"] * alpha
        b = 0.07780 * R * param_list["Tc"] / param_list["Pc"]

        # 混合规则
        a_mix = np.sum(np.outer(mole_frac, mole_frac) * np.sqrt(np.outer(a, a)))
        b_mix = np.sum(mole_frac * b)

        # 定义Peng-Robinson方程
        A = a_mix * pressure / (R * temperature) ** 2
        B = b_mix * pressure / (R * temperature)

        coeff = [1,
                 B - 1,
                 A - 3 * B ** 2 - 2 * B,
                 B ** 3 + B ** 2 - A * B]

        Z_roots = np.roots(coeff)

        Z_real = Z_roots[np.isreal(Z_roots)].real
        Z_real = np.sort(Z_real)

        densities = pressure / (Z_real * R * temperature)

        return densities, Z_real

    def calculate_molar_density_pure(self, temperature, pressure, params):
        pass

    def calculate_polymer_density(self, temperature):
        pass


class PengRobinsonMethod(PropertiesMethod):

    def calculate_molar_density_mixture(self, mole_frac, temperature, pressure):
        pass
