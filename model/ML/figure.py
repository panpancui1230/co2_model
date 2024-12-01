import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 加载数据
file_path = './0_top100_CO2.csv'
data = pd.read_csv(file_path)

# 创建 CO2 梯度分类
bins = np.linspace(data['CO2'].min(), data['CO2'].max(), 4)
labels = ['Low', 'Medium', 'High']
data['CO2_category'] = pd.cut(data['CO2'], bins=bins, labels=labels, include_lowest=True)

# 定义配色：Low -> Green, Medium -> Orange, High -> Blue
palette = {'Low': '#2ca02c',  # 绿色
           'Medium': '#ff7f0e',  # 橙色
           'High': '#1f77b4'}  # 蓝色

# 使用 pairplot 绘制图形
sns.pairplot(data=data, 
             vars=['ratio_absorb', 'PSII_content_per_leaf', 'PSII_antenna_size', 'P700_red_initial'], 
             hue='CO2_category', 
             palette=palette,
             diag_kind='kde', 
             plot_kws={'alpha': 0.7, 's': 40})

# 添加标题
plt.suptitle('Impact of Variables on CO2 Levels', y=1.02, fontsize=16)

# 显示图形
plt.show()