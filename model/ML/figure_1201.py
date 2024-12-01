import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np

# 加载数据
file_path = './500_CO2.csv'
data = pd.read_csv(file_path)
data["P700_red_initial"] *= 0.7

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.size"] = 22

# 筛选 PSII_content_per_leaf = 0.84 的数据
filtered_data = data[data['PSII_content_per_leaf'] == 3.15]

# Normalize CO2 values for color mapping based on全数据范围
norm = mcolors.Normalize(vmin=data['CO2'].min(), vmax=data['CO2'].max())
cmap = plt.cm.coolwarm

# 获取 unique 的 ratio_absorb 值，并确保有 10 个值
unique_ratios = filtered_data['ratio_absorb'].unique()
if len(unique_ratios) > 10:
    unique_ratios = unique_ratios[:10]

# 创建图表，3 行 4 列子图网格
fig, axes = plt.subplots(3, 4, figsize=(18, 12))
axes = axes.flatten()

# 绘制每个子图
for i, ratio in enumerate(unique_ratios):
    ax = axes[i]
    subset = filtered_data[filtered_data['ratio_absorb'] == ratio]
    
    scatter = ax.scatter(
        subset['PSII_antenna_size'], 
        subset['P700_red_initial'], 
        c=subset['CO2'],  # 使用 CO2 实际值作为颜色映射
        cmap=cmap, 
        norm=norm, 
        alpha=0.7, 
        s=50
    )
    
    ax.set_title(r'$Abs_{ratio}$ = ' + f'{ratio:.2f}', fontsize=20)
    # 移除每个小图的 x 和 y 标签
    ax.set_xlabel('')
    ax.set_ylabel('')

# 移除多余的子图（空白位置）
for j in range(len(unique_ratios), len(axes)):
    fig.delaxes(axes[j])

# 设置全局的 x 和 y 标签
# fig.supxlabel('$PSII_{antenna}$', fontsize=25)
# fig.supylabel('PSI/leaf area', fontsize=25)
# 设置全局的 x 和 y 标签，并调整位置
fig.supxlabel('$PSII_{antenna}$', fontsize=24, x=0.5, y=0.02)  # 调整 x 和 y 参数
fig.supylabel('PSI/leaf area', fontsize=24, x=0.01, y=0.5)  # 调整 x 和 y 参数

# 添加全局颜色条在第三行的第三列
colorbar_position = [0.60, 0.2, 0.3, 0.03]  # [x, y, width, height]
cbar_ax = fig.add_axes(colorbar_position)  # 调节位置和大小
cbar = fig.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=cmap), 
    cax=cbar_ax, 
    orientation='horizontal'  # 设置为横向
)

# 映射颜色条刻度为 CO2 实际值
co2_ticks = np.linspace(data['CO2'].min(), data['CO2'].max(), num=5)  # 对应 CO2 的实际值范围
cbar.set_ticks(co2_ticks)  # 设置颜色条刻度位置
cbar.set_ticklabels([f'{tick:.2f}' for tick in co2_ticks])  # 格式化标签为实际值
cbar.set_label('A (µmol $m^{-2} s^{-1}$)', fontsize=24)

# 调整布局
plt.tight_layout(rect=[0, 0, 1, 0.95])  # 为颜色条和子图留出空间
plt.show()
fig.savefig('1201_3.15.jpg')