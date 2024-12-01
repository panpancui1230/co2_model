import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np

# 加载数据
file_path = './500_CO2.csv'
data = pd.read_csv(file_path)
data["P700_red_initial"] *= 0.7

# 筛选 PSII_content_per_leaf = 1.61 的数据
filtered_data = data[data['PSII_content_per_leaf'] == 1.61]

# 获取 unique 的 ratio_absorb 值，并确保有 10 个值
unique_ratios = filtered_data['ratio_absorb'].unique()
if len(unique_ratios) > 10:
    unique_ratios = unique_ratios[:10]

# 计算 CO2 的百分位排名
filtered_data['CO2_percentile'] = filtered_data['CO2'].rank(pct=True)

# Normalize CO2 values for color mapping
norm = mcolors.Normalize(vmin=0, vmax=1)
cmap = plt.cm.coolwarm

# 创建图表，3 行 4 列子图网格（最后一列可能空白）
fig, axes = plt.subplots(3, 4, figsize=(18, 12))  # 3 行 4 列网格布局
axes = axes.flatten()

# 绘制每个子图
for i, ratio in enumerate(unique_ratios):
    ax = axes[i]
    subset = filtered_data[filtered_data['ratio_absorb'] == ratio]
    
    scatter = ax.scatter(
        subset['PSII_antenna_size'], 
        subset['P700_red_initial'], 
        c=subset['CO2_percentile'],  # 使用百分位排名作为颜色映射
        cmap=cmap, 
        norm=norm, 
        alpha=0.7, 
        s=50
    )
    
    ax.set_title(f'Ratio Absorb = {ratio:.2f}', fontsize=12)
    ax.set_xlabel('PSII Antenna Size', fontsize=10)
    ax.set_ylabel('P700 Red Initial', fontsize=10)

# 移除多余的子图（空白位置）
for j in range(len(unique_ratios), len(axes)):
    fig.delaxes(axes[j])

# 添加全局颜色条
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # 颜色条位置和大小 [x, y, width, height]
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap), 
    cax=cbar_ax, 
    orientation='vertical'
)

# 映射颜色条刻度为 CO2 实际值
co2_ticks = np.linspace(0, 1, num=5)  # 百分位范围
co2_actual_ticks = np.percentile(filtered_data['CO2'], co2_ticks * 100)  # 对应的实际 CO2 值
cbar.set_ticks(co2_ticks)  # 设置颜色条刻度位置
cbar.set_ticklabels([f'{tick:.2f}' for tick in co2_actual_ticks])  # 将刻度标签设置为实际值
cbar.set_label('CO2 (µmol $m^{-2} s^{-1}$)', fontsize=14)

# 调整布局
plt.tight_layout(rect=[0, 0, 0.9, 1])  # 为颜色条和子图留出空间
plt.show()