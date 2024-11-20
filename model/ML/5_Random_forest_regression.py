import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

# 读取数据
file_path = './100_CO2.csv'
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 数据集划分 (80% 训练集, 20% 测试集)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 构建随机森林回归模型
rf = RandomForestRegressor(n_estimators=100, random_state=42)

# 模型训练
rf.fit(X_train, y_train)

# 模型预测
y_pred = rf.predict(X_test)

# 模型评估
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error (MSE): {mse:.4f}")
print(f"R-squared (R2): {r2:.4f}")

# 可选：输出特征重要性
feature_importances = pd.DataFrame({
    'Feature': X.columns,
    'Importance': rf.feature_importances_
}).sort_values(by='Importance', ascending=False)

print("\nFeature Importances:")
print(feature_importances)


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