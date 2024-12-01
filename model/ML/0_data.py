#取CO2最大的10组
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from seaborn import load_dataset
file_path = './500_CO2.csv' 
# file_path='./0_top100_CO2.csv' 
df = pd.read_csv(file_path)
df["P700_red_initial"] *= 0.7
# df["PSI/PSII"] = df["P700_red_initial"] / df["PSII_content_per_leaf"] 
# df["PSII/PSI"] =  df["PSII_content_per_leaf"] /df["P700_red_initial"]
# plt.rcParams["figure.dpi"] = 300
plt.rcParams['font.sans-serif'] ='Arial'
plt.rcParams["font.size"] = 18
# plt.rcParams['mathtext.default'] = 'it'

# # output_dir = "/Users/panpan/Desktop/co2_model/model/ML/plot"
# #
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 找到 CO2 较大的 100 组值
top_20 = df.nlargest(200, "CO2")
top_20_sorted = top_20.sort_values(by="CO2",ascending=False)


output_path = './0_top200_CO2.csv'
top_20_sorted.to_csv(output_path,index=False)

print(f"文件保存至:{output_path}")


# 可视化变量对 CO2 的影响
# for col in X.columns:
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(x=df[col], y=y, edgecolor="w", s=70,color="green")
#     plt.title(f"Effect of {col} on CO2", fontsize=16)
#     plt.xlabel(col, fontsize=14)
#     plt.ylabel("CO2", fontsize=14)
#     plt.grid(False)
#     plt.show()

# # 单变量关系_组图
# # fig, axes = plt.subplots(2, 2, figsize=(12, 8))
# fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# cols = X.columns
# axes = axes.ravel()  


# # 每个特征的散点图
# custom_x_labels = [
#     '$Abs_{ratio}$', 
#     'PSII/leaf area (µmol $m^{-2} s^{-1}$)', 
#     # 'PSII $antenna$', 
#     '$PSII_{antenna}$',
#     'PSI/leaf area (µmol $m^{-2} s^{-1}$)'
# ]

# custom_xticks = [
#     [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],  # 图1
#     [0, 1, 2, 3, 4, 5, 6, 7],                      # 图2
#     [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],  # 图3
#     [0, 1, 2, 3, 4, 5, 6, 7]                       # 图4
# ]
# # for i, col in enumerate(cols):
# for i, (col, label,xticks) in enumerate(zip(cols, custom_x_labels,custom_xticks)):
#     sns.scatterplot(
#         x=df[col],
#         y=y,
#         edgecolor="w",
#         s=70,
#         color="green",
#         ax=axes[i] 
#     )
#     # axes[i].set_title(f"Effect of {col} on CO2", fontsize=14)  
#     axes[i].set_xlabel(label, fontsize=18)  
#     axes[i].set_ylabel('A (µmol $m^{-2} s^{-1}$)', fontsize=18)
#     axes[i].set_xticks(xticks)  # 设置自定义x轴刻度 
#     axes[i].grid(False)  
    
#     # if col == "P700_red_initial":
#     #     # 修改x轴刻度值
#     #     original_xticks = axes[i].get_xticks()
#     #     new_xticks = [tick * 0.7 for tick in original_xticks]
#     #     axes[i].set_xticks(new_xticks)
# plt.tight_layout()
# plt.show()

# #3. 4个变量与CO2的相关性
# correlation_matrix = df.corr()

# variables = ["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]
# co2_correlations = correlation_matrix["CO2"].loc[variables]

# # 自定义x轴变量名称
# custom_x_labels = [
#     '$Abs_{ratio}$', 
#     'PSII/leaf area', 
#     '$PSII_{antenna}$',
#     'PSI/leaf area',
# ]

# plt.figure(figsize=(10, 6))

# single_color = sns.color_palette("coolwarm", n_colors=5)[0]
# bar_plot = sns.barplot(x=custom_x_labels, y=co2_correlations.values, color=single_color)
# # bar_plot = sns.barplot(x=custom_x_labels, y=co2_correlations.values, palette="coolwarm")
# # plt.title("Impact of Variables on CO2", fontsize=16)
# # plt.xlabel("Variables", fontsize=18)
# plt.ylabel("Correlation with A", fontsize=18)
# plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  
# # plt.grid(axis="y", linestyle="--", alpha=0.7)
# plt.grid(False)

# # 调整 y 轴范围，增高边框线
# y_min, y_max = plt.ylim()
# plt.ylim(y_min-0.1, y_max + 0.1) 
# # 坐标刻度字体大小
# plt.xticks(fontsize=18)  # 减小x轴刻度字体大小
# plt.yticks(fontsize=18)  # 减小y轴刻度字体大小
# plt.xlabel(None)

# ax = plt.gca()  # 获取当前轴
# ax.tick_params(axis='x', pad=14)  # 增加 x 轴刻度字体与轴的距离


# for index, value in enumerate(co2_correlations.values):
#     plt.text(
#         x=index, 
#         y=value + 0.02 if value > 0 else value - 0.07,  # 标签位置稍高于柱形顶部
#         s=f"{value:.2f}", 
#         ha='center', va='bottom', fontsize=18, color="black", 
#         bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')  
#     )

# plt.show()


# #top 100
# file_path = './0_top100_CO2.csv' 
# df = pd.read_csv(file_path)

# sns.pairplot(df.iloc[:, :4], diag_kind='kde', corner=True, plot_kws={'alpha': 0.9})

# plt.suptitle("Scatter Plot Matrix of Variables", fontsize=16, y=1.02)
# plt.show()

# # 加载企鹅数据集
# penguins = sns.load_dataset("penguins")

# # 重命名 penguins 数据集的列，映射到你提供的列名
# penguins = penguins.rename(columns={
#     "bill_length_mm": "ratio_absorb",
#     "bill_depth_mm": "PSII_content_per_leaf",
#     "flipper_length_mm": "PSII_antenna_size",
#     "body_mass_g": "P700_red_initial"
# })

# # 删除缺失值
# penguins = penguins.dropna(subset=["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"])

# # 绘制 Pairplot
# g=sns.pairplot(
#     penguins[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial", "species"]],
#     hue="species",
#     diag_kind="kde",
#     corner=True,
#     palette="Set2",
#     plot_kws={'alpha': 0.6}
# )

# # 动态设置每个轴的刻度值
# # 动态调整刻度
# for ax in g.axes.flat:
#     if ax is not None:
#         xlim = ax.get_xlim()
#         ylim = ax.get_ylim()
#         ax.set_xticks(np.linspace(xlim[0], xlim[1], 5))  # 设置5个均匀分布的刻度
#         ax.set_yticks(np.linspace(ylim[0], ylim[1], 5))

# # 添加标题
# plt.suptitle("Pairplot of Penguins Dataset (Grouped by Species)", y=1.02, fontsize=16)

# # 显示图像
# # plt.show()
# plt.show()
# # print(f"总点数: {df.shape[0]}")

# # # 绘制每个变量的分布
# # for i, col in enumerate(df.columns[:4], 1):  # 遍历前4列（自变量）
# #     plt.subplot(2, 2, i)  # 2行2列的子图布局
# #     sns.histplot(df[col], kde=True, color='blue', bins=20)  # 直方图 + KDE
# #     plt.title(f"Distribution of {col}")
# #     plt.xlabel(col)
# #     plt.ylabel("Frequency")

# # plt.tight_layout()
# # plt.show()

