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


