import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

# 读取数据
file_path = './100_CO2.csv'
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 数据集划分 (80% 训练集, 20% 测试集)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 构建随机森林回归模型
rf = RandomForestRegressor(n_estimators=100,random_state=42)

# 模型训练
rf.fit(X_train, y_train)

# 模型预测
y_train_pred = rf.predict(X_train)
y_test_pred = rf.predict(X_test)

train_mse = mean_squared_error(y_train, y_train_pred)  # 训练集 MSE
train_rmse = train_mse ** 0.5                          # 训练集 RMSE
train_mae = mean_absolute_error(y_train, y_train_pred) # 训练集 MAE
train_r2 = r2_score(y_train, y_train_pred)             # 训练集 R2

# 测试集评估
test_mse = mean_squared_error(y_test, y_test_pred)     # 测试集 MSE
test_rmse = test_mse ** 0.5                            # 测试集 RMSE
test_mae = mean_absolute_error(y_test, y_test_pred)    # 测试集 MAE
test_r2 = r2_score(y_test, y_test_pred)                # 测试集 R2

# 输出训练集评估结果
print("\nTraining Data Evaluation:")
print(f"Train Mean Squared Error (MSE): {train_mse:.4f}")
print(f"Train Root Mean Squared Error (RMSE): {train_rmse:.4f}")
print(f"Train Mean Absolute Error (MAE): {train_mae:.4f}")
print(f"Train R-squared (R2): {train_r2:.4f}")

# 输出测试集评估结果
print("\nTesting Data Evaluation:")
print(f"Test Mean Squared Error (MSE): {test_mse:.4f}")
print(f"Test Root Mean Squared Error (RMSE): {test_rmse:.4f}")
print(f"Test Mean Absolute Error (MAE): {test_mae:.4f}")
print(f"Test R-squared (R2): {test_r2:.4f}")



#残差图 (Residual Plot)
# Calculate residuals
residuals = y_test - y_test_pred

# Create a DataFrame for seaborn
residuals_df = pd.DataFrame({
    "Predicted Values": y_test_pred,
    "Residuals": residuals
})

# Plot residuals using seaborn
plt.figure(figsize=(8, 6))
sns.scatterplot(data=residuals_df, x="Predicted Values", y="Residuals", color='green', alpha=0.7, label="Residuals")
plt.axhline(y=0, color='r', linestyle='--', label="Residual = 0")  # Horizontal line at y=0

# Add labels and title
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.legend(loc='upper left')
plt.grid(False)
plt.show()


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


import matplotlib.pyplot as plt

# 定义评估指标和对应的值
metrics = ["MSE", "RMSE", "MAE", "R^2"]
train_values = [train_mse, train_rmse, train_mae, train_r2]
test_values = [test_mse, test_rmse, test_mae, test_r2]

# 创建折线图
plt.figure(figsize=(10, 6))
plt.plot(metrics, train_values, marker='o', label="Training Data", linestyle='-', linewidth=2, markersize=8, color='blue')
plt.plot(metrics, test_values, marker='o', label="Testing Data", linestyle='-', linewidth=2, markersize=8, color='green')

# 添加数值标签
for i, value in enumerate(train_values):
    plt.text(i, value + 0.01, f"{value:.4f}", ha="center", va="bottom", fontsize=10, color='blue')
for i, value in enumerate(test_values):
    plt.text(i, value + 0.01, f"{value:.4f}", ha="center", va="bottom", fontsize=10, color='green')

# 设置标题和轴标签
plt.title("Model Evaluation Metrics: Training vs Testing", fontsize=16)
plt.ylabel("Metric Value", fontsize=14)
plt.xlabel("Metrics", fontsize=14)

# 添加网格和图例
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.legend(fontsize=12)

# 显示图表
plt.tight_layout()
plt.show()


