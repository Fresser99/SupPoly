from abc import ABC, abstractmethod


class PropertiesMethod(ABC):

    @abstractmethod
    def calculate_molar_density_mixture(self, *args, **kwargs):
        pass

    @abstractmethod
    def calculate_molar_density_pure(self, *args, **kwargs):
        pass


class UserMethod(PropertiesMethod):

    def calculate_molar_density_mixture(self, mole_frac, temperature, pressure, param_list):
        return 0.

    def calculate_molar_density_pure(self, temperature, pressure, params):
        return 10.

class PengRobinsonMethod(PropertiesMethod):

    def calculate_molar_density_mixture(self, mole_frac,temperature,pressure):

        return 10.