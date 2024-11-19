import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# 从 CSV 文件中读取数据
file_path = './100_CO2.csv'  # 请将此路径替换为实际 CSV 文件的路径
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 决策树回归模型训练
regr = DecisionTreeRegressor(max_depth=10)  # 设置最大深度为 2
# regr = DecisionTreeRegressor(
#     max_depth=5,               # 增加最大深度
#     min_samples_split=10,      # 设置分裂所需的最少样本数
#     min_samples_leaf=5         # 每个叶子节点的最少样本数
# )
regr.fit(X, y)


# 测试数据（从训练数据中提取一部分或生成新的数据）
X_test = pd.DataFrame({
    "ratio_absorb": [0.84],
    "PSII_content_per_leaf": [0.7],
    "PSII_antenna_size": [0.5],
    "P700_red_initial": [0.7],
})

# 进行预测
predictions = regr.predict(X_test)
print("Predicted CO2 values:", predictions)

# 评估模型性能
# 使用训练数据进行预测
y_train_pred = regr.predict(X)

# 计算评估指标
mse = mean_squared_error(y, y_train_pred)
rmse = mse ** 0.5
mae = mean_absolute_error(y, y_train_pred)
r2 = r2_score(y, y_train_pred)

# 输出评估结果
print("\nModel Evaluation on Training Data:")
print(f"Mean Squared Error (MSE): {mse:.4f}")
print(f"Root Mean Squared Error (RMSE): {rmse:.4f}")
print(f"Mean Absolute Error (MAE): {mae:.4f}")
print(f"R^2 Score: {r2:.4f}")