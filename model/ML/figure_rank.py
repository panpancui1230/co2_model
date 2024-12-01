import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd

# 加载数据
file_path = './500_CO2.csv'
data = pd.read_csv(file_path)
data["P700_red_initial"] *= 0.7
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.size"] = 18

# 将 CO2 转换为百分位排名
data['CO2_percentile'] = data['CO2'].rank(pct=True)

# Normalize CO2 percentile values for color mapping
norm = mcolors.Normalize(vmin=0, vmax=1)  # 因为是百分位，范围固定为 0 到 1
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

# Map scatter plots with color based on CO2 percentile
g.map_offdiag(lambda x, y, **kwargs: plt.scatter(x, y, c=data['CO2_percentile'], cmap=cmap, norm=norm, alpha=0.7, s=40))

# Map histograms to the diagonal
# g.map_diag(sns.histplot, kde=False, color='#f0c76e', edgecolor=None)
g.map_diag(sns.histplot, kde=False, color='#faba4b', edgecolor=None)
# g.map_diag(sns.histplot, kde=False, color='#ebdb2f', edgecolor=None)
#fabe55

# Add colorbar to represent CO2 values
cbar = g.fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=g.axes, orientation='vertical', shrink=0.8)
cbar.set_label('A (µmol $m^{-2} s^{-1}$)')

# 设置色条的刻度为百分位对应的实际 CO2 值
percentile_ticks = [0.0, 0.25, 0.5, 0.75, 1.0]  # 百分位位置
actual_co2_values = np.percentile(data['CO2'], [p * 100 for p in percentile_ticks])  # 转换为实际 CO2 值

cbar.set_ticks(percentile_ticks)  # 设置色条的刻度位置
cbar.set_ticklabels([f'{val:.2f}' for val in actual_co2_values])  # 对应实际 CO2 值作为标签

# 添加自定义标签并减小字体
for i, label in enumerate(custom_x_labels):
    g.axes[i, 0].set_ylabel(label, fontsize=18, labelpad=10)  # 设置 y 轴标签字体大小
    g.axes[-1, i].set_xlabel(label, fontsize=18, labelpad=10)  # 设置 x 轴标签字体大小

# 设置 x 和 y 轴刻度（保留一位小数）
tick_settings = {
    0: [0.2, 0.4, 0.6, 0.8],  # 第 1 个变量
    1: [1.0, 3.0, 5.0, 7.0],  # 第 2 个变量
    2: [0.2, 0.4, 0.6, 0.8],  # 第 3 个变量
    3: [1.0, 3.0, 5.0, 7.0]   # 第 4 个变量
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
g.savefig('rank.jpg')