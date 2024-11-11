import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm


class PolymerMWD:
    def __init__(self, M0, M1, M2):
        """
        初始化聚合物分子量分布重构器

        参数:
        M0: 零阶矩 (总摩尔数)
        M1: 一阶矩
        M2: 二阶矩
        """
        self.M0 = M0
        self.M1 = M1
        self.M2 = M2

        # 计算基本参数
        self.Mn = self.M1 / self.M0  # 数均分子量
        self.Mw = self.M2 / self.M1  # 重均分子量
        self.PDI = self.Mw / self.Mn  # 多分散性指数
        # 计算对数正态分布参数
        self.sigma = np.sqrt(np.log(self.PDI))
        self.mu = np.log(self.Mn) - (self.sigma ** 2) / 2

    def calculate_distribution(self, M):
        """
        计算给定分子量点的分布值

        参数:
        M: 分子量值或数组

        返回:
        w: 对应的权重分数
        """
        w = (1 / (M * self.sigma * np.sqrt(2 * np.pi))) * \
            np.exp(-(np.log(M) - self.mu) ** 2 / (2 * self.sigma ** 2))
        return w

    def plot_distribution(self, M_range=None, num_points=1000):
        """
        绘制分子量分布

        参数:
        M_range: 分子量范围元组 (min_M, max_M)
        num_points: 计算点数
        """
        if M_range is None:
            # 默认范围：Mn/10 到 Mw*5
            M_range = (1, self.Mw * 5)

            # 生成对数均匀分布的分子量点
        M = np.logspace(np.log10(M_range[0]), np.log10(M_range[1]), num_points)
        w = self.calculate_distribution(M)

        # 创建图形
        plt.figure(figsize=(10, 6))
        # 绘制分布曲线
        plt.semilogx(M, w, 'b-', label='MWD')

        # 添加垂直线标记 Mn 和 Mw
        plt.axvline(x=self.Mn, color='r', linestyle='--', label=f'Mn = {self.Mn:.0f}')
        plt.axvline(x=self.Mw, color='g', linestyle='--', label=f'Mw = {self.Mw:.0f}')

        # 设置图形属性
        plt.xlabel('Mn (g/mol)')
        plt.ylabel('Weight')
        plt.title(f'MWD (PDI = {self.PDI:.2f})')
        plt.grid(True, which="both", ls="-", alpha=0.2)
        plt.legend()

        # 添加分布参数文本
        text = f'M₀ = {self.M0:.2e} mol\n' \
               f'M₁ = {self.M1:.2e} g/mol\n' \
               f'M₂ = {self.M2:.2e} (g/mol)²\n' \
               f'Mn = {self.Mn:.2e} g/mol\n' \
               f'Mw = {self.Mw:.2e} g/mol\n' \
               f'PDI = {self.PDI:.2f}'
        plt.text(0.02, 0.98, text, transform=plt.gca().transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        plt.show()

    def get_distribution_parameters(self):
        """
        返回分布的关键参数
        """
        return {
            'M0': self.M0,
            'M1': self.M1,
            'M2': self.M2,
            'Mn': self.Mn,
            'Mw': self.Mw,
            'PDI': self.PDI,
            'mu': self.mu,
            'sigma': self.sigma
        }
