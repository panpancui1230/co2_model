# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# # 读取上传的文件
# file_path = '/Users/panpan/Desktop/co2_model/model/ML/plot/0_top20_CO2.csv'
# df = pd.read_csv(file_path)

# # 确保所有列为数值类型，并删除非数值或缺失值
# df = df.apply(pd.to_numeric, errors='coerce').dropna()

# # 提取前4个变量
# variables = df.iloc[:, :4]

# # 在一张图上绘制所有变量的分布
# plt.figure(figsize=(10, 6))
# for col in variables.columns:
#     sns.kdeplot(variables[col], label=col, fill=True, alpha=0.5)

# # 添加标题和图例
# plt.title("Distributions of Variables", fontsize=16)
# plt.xlabel("Value", fontsize=14)
# plt.ylabel("Density", fontsize=14)
# plt.legend(title="Variables", fontsize=12)
# plt.tight_layout()
# plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取上传的文件
file_path = '/Users/panpan/Desktop/co2_model/model/ML/plot/0_top100_CO2.csv'
df = pd.read_csv(file_path)

# 为数据添加标记列，用于区分变量
df_melted = df.iloc[:, :4].melt(var_name="Variable", value_name="Value")

# 绘制散点图矩阵
sns.pairplot(
    df.iloc[:, :4], 
    diag_kind="kde", 
    corner=True, 
    hue="Variable",  # 通过变量区分颜色
    palette="Set2",  # 使用对比鲜明的调色板
    plot_kws={"alpha": 0.7}  # 设置透明度
)

# 设置标题
plt.suptitle("Scatter Plot Matrix with Different Colors per Variable", y=1.02, fontsize=16)
plt.show()

