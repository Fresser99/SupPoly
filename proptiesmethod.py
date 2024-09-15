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
        pass

    def calculate_molar_density_pure(self, temperature, pressure, params):
        pass

class PengRobinsonMethod(PropertiesMethod):

    def calculate_molar_density_mixture(self, mole_frac,temperature,pressure):
        pass