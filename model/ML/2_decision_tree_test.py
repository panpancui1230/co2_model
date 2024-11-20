import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# 从 CSV 文件中读取数据
file_path = './100_CO2.csv'  
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

#划分训练和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 决策树回归模型训练
# regr = DecisionTreeRegressor(max_depth=20)

regr = DecisionTreeRegressor(
    max_depth=20,               # 增加最大深度
    min_samples_split=3,      # 设置分裂所需的最少样本数
    min_samples_leaf=1         # 每个叶子节点的最少样本数
)
regr.fit(X_train, y_train)

# 使用测试数据进行预测
y_test_pred = regr.predict(X_test)

# 计算测试集上的评估指标
test_mse = mean_squared_error(y_test, y_test_pred)
test_rmse = test_mse ** 0.5
test_mae = mean_absolute_error(y_test, y_test_pred)
test_r2 = r2_score(y_test, y_test_pred)

# 输出测试集评估结果
print("Model Evaluation on Testing Data:")
print(f"Test MSE: {test_mse:.4f}")
# print(f"Test RMSE: {test_rmse:.4f}")
# print(f"Test MAE: {test_mae:.4f}")
print(f"Test R^2 Score: {test_r2:.4f}")


import matplotlib.pyplot as plt
import seaborn as sns

# 创建条形图以可视化评估指标，使用配色 palette="coolwarm"
metrics = ["MSE", "RMSE", "MAE", "R^2"]
values = [test_mse, test_rmse, test_mae, test_r2]

# 设置配色方案
palette = sns.color_palette("coolwarm", len(metrics))

# 创建图形
plt.figure(figsize=(10, 6))
bars = plt.bar(metrics, values, color=palette)

# 添加数值标签
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2, height + 0.01, f"{height:.4f}", ha="center", va="bottom", fontsize=12)

# 设置标题和轴标签
plt.title("Model Evaluation Metrics", fontsize=16)
plt.ylabel("Metric Value", fontsize=14)
plt.xlabel("Metrics", fontsize=14)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# 显示图表
plt.tight_layout()
plt.show()

# #残差图 (Residual Plot)
# # Calculate residuals
# residuals = y_test - y_test_pred

# # Create a DataFrame for seaborn
# residuals_df = pd.DataFrame({
#     "Predicted Values": y_test_pred,
#     "Residuals": residuals
# })

# # Plot residuals using seaborn
# plt.figure(figsize=(8, 6))
# sns.scatterplot(data=residuals_df, x="Predicted Values", y="Residuals", color='green', alpha=0.7, label="Residuals")
# plt.axhline(y=0, color='r', linestyle='--', label="Residual = 0")  # Horizontal line at y=0

# # Add labels and title
# plt.xlabel("Predicted Values")
# plt.ylabel("Residuals")
# plt.title("Residual Plot")
# plt.legend(loc='upper left')
# plt.grid(False)
# plt.show()


#实际值 vs 预测值散点图
# Construct a DataFrame for use with seaborn
results_df = pd.DataFrame({
    "Actual": y_test,
    "Predicted": y_test_pred
})

# Plot scatterplot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=results_df, x="Actual", y="Predicted", color='blue', alpha=0.7, label="Predicted vs Actual")
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', label="Ideal Prediction")  # Ideal prediction line

# Add legend and title
plt.xlabel("Actual Values")
plt.ylabel("Predicted Values")
plt.title("Predicted vs Actual Values")
plt.legend(loc='upper left')
plt.grid(False)
plt.show()