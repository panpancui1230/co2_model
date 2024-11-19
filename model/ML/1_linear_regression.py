import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 加载数据
file_path = './100_CO2.csv'

data = pd.read_csv(file_path)

# 自变量和因变量
X = data.iloc[:, :-1]  # 前4列为自变量
y = data.iloc[:, -1]   # 最后一列为因变量

# 训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 线性回归模型
lr_model = LinearRegression()

# 训练模型
lr_model.fit(X_train, y_train)

# 测试集预测
y_pred = lr_model.predict(X_test)

# 评估模型性能
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

# 性能评估
print(f"均方误差 (MSE): {mse}")
print(f"决定系数 (R²): {r2}")

# 模型的系数和截距
coefficients = lr_model.coef_  # 每个自变量的系数
intercept = lr_model.intercept_  # 截距

print("模型系数:", coefficients)
print("模型截距:", intercept)

# 预测新数据
new_data = [[0.84, 0.7, 0.5, 0.7]] 
predicted_value = lr_model.predict(new_data)
print(f"预测值: {predicted_value[0]}")


# #实际值 vs 预测值散点图
# # Construct a DataFrame for use with seaborn
# results_df = pd.DataFrame({
#     "Actual": y_test,
#     "Predicted": y_pred
# })

# # Plot scatterplot
# plt.figure(figsize=(8, 6))
# sns.scatterplot(data=results_df, x="Actual", y="Predicted", color='blue', alpha=0.7, label="Predicted vs Actual")
# plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', label="Ideal Prediction")  # Ideal prediction line

# # Add legend and title
# plt.xlabel("Actual Values")
# plt.ylabel("Predicted Values")
# plt.title("Predicted vs Actual Values")
# plt.legend(loc='upper left')
# plt.grid(False)
# plt.show()


# #残差图 (Residual Plot)
# # Calculate residuals
# residuals = y_test - y_pred

# # Create a DataFrame for seaborn
# residuals_df = pd.DataFrame({
#     "Predicted Values": y_pred,
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

# #模型系数的条形图
# output_dir = "/Users/panpan/Desktop/co2_model/model/ML/plot"

# # Create a DataFrame for coefficients
# coefficients_df = pd.DataFrame({
#     "Features": X.columns,  # Feature names
#     "Coefficients": coefficients  # Corresponding coefficients
# })

# # Plot coefficients using seaborn
# plt.figure(figsize=(8, 6))
# sns.barplot(data=coefficients_df, x="Features", y="Coefficients", color='green', alpha=0.7)

# # Add labels and title
# plt.xlabel("Features")
# plt.ylabel("Coefficient Values")
# plt.title("Model Coefficients Visualization")
# plt.grid(False)

# output_path = os.path.join(output_dir, "模型系数的条形图.png")
# plt.savefig(output_path, dpi=300, bbox_inches="tight")  # 高分辨率保存图表
# plt.show()

plt.figure(figsize=(8, 6))
plt.plot(range(len(y_test)), y_test, 'bo-', label='实际值')
plt.plot(range(len(y_pred)), y_pred, 'r*-', label='预测值')
plt.xlabel('样本索引')
plt.ylabel('值')
plt.title('实际值与预测值趋势对比')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()