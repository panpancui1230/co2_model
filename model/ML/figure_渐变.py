import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd

# 加载数据
file_path = './500_CO2.csv'
data = pd.read_csv(file_path)
data["P700_red_initial"] *= 0.7
# plt.rcParams["figure.dpi"] = 300
plt.rcParams['font.sans-serif'] ='Arial'
plt.rcParams["font.size"] = 18


# Normalize CO2 values for color mapping
norm = mcolors.Normalize(vmin=data['CO2'].min(), vmax=data['CO2'].max())
cmap = plt.cm.coolwarm

# 自定义标签
custom_x_labels = [
    '$Abs_{ratio}$', 
    'PSII/leaf area', 
    '$PSII_{antenna}$',
    'PSI/leaf area'
]

# Create PairGrid
g = sns.PairGrid(data, vars=['ratio_absorb', 'PSII_content_per_leaf', 'PSII_antenna_size', 'P700_red_initial'])

# Map scatter plots with color based on CO2
g.map_offdiag(lambda x, y, **kwargs: plt.scatter(x, y, c=data['CO2'], cmap=cmap, norm=norm, alpha=0.7, s=40))

# Map histograms to the diagonal
# g.map_diag(sns.histplot, kde=False, color='#f2d094', edgecolor=None)
g.map_diag(sns.histplot, kde=False, color='#f0c76e', edgecolor=None)
# Add colorbar to represent CO2 values
cbar = g.fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=g.axes, orientation='vertical', shrink=0.8)
cbar.set_label('A (µmol $m^{-2} s^{-1}$)')


# 添加自定义标签并减小字体
for i, label in enumerate(custom_x_labels):
    g.axes[i, 0].set_ylabel(label, fontsize=18, labelpad=10)  # 设置 y 轴标签字体大小
    g.axes[-1, i].set_xlabel(label, fontsize=18, labelpad=10)  # 设置 x 轴标签字体大小

# 设置 x 和 y 轴刻度（保留一位小数）
tick_settings = {
    0: [0.2, 0.4, 0.6,0.8],  # 第 1 个变量
    1: [1.0, 3.0, 5.0,7.0],  # 第 2 个变量
    2: [0.2, 0.4, 0.6,0.8],  # 第 3 个变量
    3: [1.0, 3.0, 5.0,7.0]   # 第 4 个变量
}

# 应用刻度和格式化标签
for i, ticks in tick_settings.items():
    # x 轴设置
    g.axes[-1, i].set_xticks(ticks)
    g.axes[-1, i].set_xticklabels([f'{tick:.1f}' for tick in ticks], fontsize=18)
    
    # y 轴设置
    g.axes[i, 0].set_yticks(ticks)
    g.axes[i, 0].set_yticklabels([f'{tick:.1f}' for tick in ticks], fontsize=18)

# Add a title
g.fig.suptitle('Pairplot of Variables with CO2 Gradient (Coolwarm Colormap)', y=1.02, fontsize=16)

# Show the plot
plt.show() 
g.savefig('111111.jpg')