#取CO2最大的10组
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# 读取数据
file_path = './100_CO2.csv'  # 替换为实际路径
df = pd.read_csv(file_path)

output_dir = "/Users/panpan/Desktop/co2_model/model/ML/plot"
# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 找到 CO2 较大的 10 组值
top_20 = df.nlargest(20, "CO2")
top_20_sorted = top_20.sort_values(by="CO2",ascending=False)

#保存为csv文件
output_path = './plot/0_top20_CO2'
top_20_sorted.to_csv(output_path,index=False)

print(f"文件保存至:{output_path}")


# 可视化变量对 CO2 的影响
# 1. 单变量关系
# for col in X.columns:
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(x=df[col], y=y, edgecolor="w", s=70,color="green")
#     plt.title(f"Effect of {col} on CO2", fontsize=16)
#     plt.xlabel(col, fontsize=14)
#     plt.ylabel("CO2", fontsize=14)
#     plt.grid(False)
#     plt.show()

# #2. 单变量关系_组图
# fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# # 列表化 X.columns 和子图的索引
# cols = X.columns
# axes = axes.ravel()  

# # 循环绘制每个特征的散点图
# for i, col in enumerate(cols):
#     sns.scatterplot(
#         x=df[col],
#         y=y,
#         edgecolor="w",
#         s=70,
#         color="green",
#         ax=axes[i] 
#     )
#     axes[i].set_title(f"Effect of {col} on CO2", fontsize=14)  
#     axes[i].set_xlabel(col, fontsize=12)  
#     axes[i].set_ylabel("CO2", fontsize=12)  
#     axes[i].grid(False)  

# # 调整布局
# plt.tight_layout()
# plt.show()

#3. 4个变量与CO2的相关性
# 计算相关性矩阵
correlation_matrix = df.corr()

# 提取4个变量对 CO2 的相关性
variables = ["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]
co2_correlations = correlation_matrix["CO2"].loc[variables]

# # 绘制条形图显示影响大小
# plt.figure(figsize=(10, 6))
# sns.barplot(x=co2_correlations.index, y=co2_correlations.values, palette="coolwarm")
# plt.title("Impact of Variables on CO2", fontsize=16)
# plt.xlabel("Variables", fontsize=14)
# plt.ylabel("Correlation with CO2", fontsize=14)
# plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # 添加零线
# plt.grid(axis="y", linestyle="--", alpha=0.7)
# plt.show()



# 绘制条形图并添加具体数值，增高边框线避免数值被挡住
plt.figure(figsize=(10, 6))
bar_plot = sns.barplot(x=co2_correlations.index, y=co2_correlations.values, palette="coolwarm")
plt.title("Impact of Variables on CO2", fontsize=16)
plt.xlabel("Variables", fontsize=14)
plt.ylabel("Correlation with CO2", fontsize=14)
plt.axhline(0, color="black", linewidth=0.8, linestyle="--")  # 添加零线
plt.grid(axis="y", linestyle="--", alpha=0.7)

# 调整 y 轴范围，增高边框线
y_min, y_max = plt.ylim()
plt.ylim(y_min-0.1, y_max + 0.1)  # 在顶部增加 0.1 的空白以容纳标签

# 添加数值标签
for index, value in enumerate(co2_correlations.values):
    plt.text(
        x=index, 
        y=value + 0.02 if value > 0 else value - 0.07,  # 标签位置稍高于柱形顶部
        s=f"{value:.2f}", 
        ha='center', va='bottom', fontsize=12, color="black", 
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')  # 添加背景框避免挡住数值
    )

plt.show()

