import pandas as pd
import matplotlib.pyplot as plt

# 数据加载
data = pd.DataFrame({
    'ratio_absorb': [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.81111111, 0.9, 0.81111111, 0.9, 0.9],
    'PSII_content_per_leaf': [1.61, 3.15, 2.38, 1.61, 2.38, 3.15, 2.38, 1.61, 1.61, 2.38, 1.61, 1.61, 3.92, 3.92, 3.92, 3.15, 0.84, 2.38, 2.38, 0.84],
    'PSII_antenna_size': [0.9, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.81111111, 0.72222222, 0.72222222, 0.72222222, 0.81111111, 0.9, 0.81111111, 0.72222222, 0.9],
    'P700_red_initial': [10.0, 10.0, 7.8, 5.6, 8.9, 8.9, 10.0, 6.7, 7.8, 6.7, 10.0, 8.9, 10.0, 8.9, 7.8, 10.0, 10.0, 7.8, 4.5, 8.9],
    'CO2': [12.73080038, 12.66728341, 12.63104185, 12.45653595, 12.32028771, 12.2284811, 12.19090454, 12.06655971, 11.75171894, 11.72904504, 11.72232874, 11.720607, 11.62710929, 11.62049643, 11.61910348, 11.61724222, 11.61457953, 11.58805619, 11.5701609, 11.55406648]
})

# 创建散点图
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    data['ratio_absorb'],  # 横轴：第一个自变量
    data['PSII_content_per_leaf'],  # 纵轴：第二个自变量
    c=data['CO2'],  # 颜色：因变量
    s=data['PSII_antenna_size'] * 200,  # 点大小：第三个自变量
    alpha=0.8,
    cmap='viridis'
)

# 添加颜色条
plt.colorbar(scatter, label='CO2')

# 设置图表标签和标题
plt.xlabel('ratio_absorb', fontsize=14)
plt.ylabel('PSII_content_per_leaf', fontsize=14)
plt.title('CO2 by Four Variables', fontsize=16)

# 添加注释
plt.tight_layout()
plt.show()