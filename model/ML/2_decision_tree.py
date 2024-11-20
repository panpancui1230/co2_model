import pandas as pd
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
regr = DecisionTreeRegressor(max_depth=10)
regr.fit(X_train, y_train)

# 使用测试数据进行预测
y_test_pred = regr.predict(X_test)

# 使用训练数据进行预测
y_train_pred = regr.predict(X_train)

# # 测试数据（从训练数据中提取一部分或生成新的数据）
# X_test = pd.DataFrame({
#     "ratio_absorb": [0.84],
#     "PSII_content_per_leaf": [0.7],
#     "PSII_antenna_size": [0.5],
#     "P700_red_initial": [0.7],
# })
# # 进行预测
# predictions = regr.predict(X_test)
# print("Predicted CO2 values:", predictions)

# 计算训练集上的评估指标
train_mse = mean_squared_error(y_train, y_train_pred)
train_rmse = train_mse ** 0.5
train_mae = mean_absolute_error(y_train, y_train_pred)
train_r2 = r2_score(y_train, y_train_pred)

# 计算测试集上的评估指标
test_mse = mean_squared_error(y_test, y_test_pred)
test_rmse = test_mse ** 0.5
test_mae = mean_absolute_error(y_test, y_test_pred)
test_r2 = r2_score(y_test, y_test_pred)

# 输出训练集评估结果
print("Model Evaluation on Training Data:")
print(f"Train MSE: {train_mse:.4f}")
print(f"Train RMSE: {train_rmse:.4f}")
print(f"Train MAE: {train_mae:.4f}")
print(f"Train R^2 Score: {train_r2:.4f}")

# 输出测试集评估结果
print("Model Evaluation on Testing Data:")
print(f"Test MSE: {test_mse:.4f}")
print(f"Test RMSE: {test_rmse:.4f}")
print(f"Test MAE: {test_mae:.4f}")
print(f"Test R^2 Score: {test_r2:.4f}")