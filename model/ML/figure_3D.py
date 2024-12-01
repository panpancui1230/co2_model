import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors

# 加载数据
file_path = './500_CO2.csv'  # 替换为实际文件路径
data = pd.read_csv(file_path)

# 数据预处理
data["PSI/PSII"] = data["P700_red_initial"]
data["P700_red_initial"] *= 0.7  # 调整 P700_red_initial
# data["PSI/PSII"] = data["P700_red_initial"] / data["PSII_content_per_leaf"]  # 计算 PSI/PSII

# 归一化 CO2 值以用于颜色映射
norm = mcolors.Normalize(vmin=data['CO2'].min(), vmax=data['CO2'].max())
cmap = plt.cm.coolwarm

# 创建 3D 散点图
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 绘制散点图，颜色基于 CO2 值
sc = ax.scatter(
    data['ratio_absorb'], 
    data['PSI/PSII'], 
    data['PSII_antenna_size'], 
    c=data['CO2'], 
    cmap=cmap, 
    norm=norm, 
    alpha=0.8, 
    s=50
)

# 设置坐标轴标签
ax.set_xlabel('Ratio Absorb', fontsize=12)
ax.set_ylabel('PSI/PSII', fontsize=12)
ax.set_zlabel('PSII Antenna Size', fontsize=12)
ax.set_title('3D Scatter Plot Colored by CO2', fontsize=14)

# 添加颜色条以表示 CO2 值
cbar = plt.colorbar(sc, shrink=0.6, aspect=10)
cbar.set_label('CO2 (µmol $m^{-2} s^{-1}$)', fontsize=12)

# 显示图形
plt.show()