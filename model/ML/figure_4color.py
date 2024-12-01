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

# 使用 Percentile 将 CO2 值分组
percentiles = np.percentile(data['CO2'], [25, 50, 75])  # 分成 4 个区间
data['CO2_group'] = pd.cut(
    data['CO2'],
    bins=[-np.inf] + list(percentiles) + [np.inf],
    labels=range(4)
)

# 从 coolwarm 色图中均匀提取 4 个颜色
cmap = plt.cm.coolwarm
color_palette = [cmap(i / 3) for i in range(4)]  # 提取 4 个等间距的颜色

# 自定义标签
custom_x_labels = [
    '$Abs_{ratio}$',
    'PSII/leaf area',
    '$PSII_{antenna}$',
    'PSI/leaf area'
]

# 创建 PairGrid
g = sns.PairGrid(data, vars=['ratio_absorb', 'PSII_content_per_leaf', 'PSII_antenna_size', 'P700_red_initial'])

# 映射散点图，按 CO2 分组上色
g.map_offdiag(lambda x, y, **kwargs: plt.scatter(
    x, y, c=data['CO2_group'].astype(int).map(lambda idx: color_palette[idx]), alpha=0.7, s=40
))

# 修改对角线柱状图的颜色
g.map_diag(sns.histplot, kde=False, color='#f2b957', edgecolor=None)

# 添加颜色图例
cbar = g.fig.colorbar(
    plt.cm.ScalarMappable(
        norm=mcolors.BoundaryNorm([0, 1, 2, 3, 4], len(color_palette)), 
        cmap=mcolors.ListedColormap(color_palette)
    ),
    ax=g.axes,
    orientation='vertical',
    shrink=0.8
)
cbar.set_label('CO2 Percentile Groups')

# 添加自定义标签并减小字体
for i, label in enumerate(custom_x_labels):
    g.axes[i, 0].set_ylabel(label, fontsize=18, labelpad=10)
    g.axes[-1, i].set_xlabel(label, fontsize=18, labelpad=10)

# 设置 x 和 y 轴刻度（保留一位小数）
tick_settings = {
    0: [0.2, 0.4, 0.6, 0.8],
    1: [1.0, 3.0, 5.0, 7.0],
    2: [0.2, 0.4, 0.6, 0.8],
    3: [1.0, 3.0, 5.0, 7.0]
}
for i, ticks in tick_settings.items():
    g.axes[-1, i].set_xticks(ticks)
    g.axes[-1, i].set_xticklabels([f'{tick:.1f}' for tick in ticks], fontsize=18)
    g.axes[i, 0].set_yticks(ticks)
    g.axes[i, 0].set_yticklabels([f'{tick:.1f}' for tick in ticks], fontsize=18)

# 添加标题
g.fig.suptitle('Pairplot of Variables with CO2 Ranking (Percentile Groups)', y=1.02, fontsize=16)

# 显示图像
plt.show()