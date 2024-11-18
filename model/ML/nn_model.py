import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error,r2_score

# 加载数据
file_path = '/Users/panpan/Desktop/co2_model/model/ML/result_100.csv'
data = pd.read_csv(file_path)

# 自变量和因变量
X = data.iloc[:,:-1]
y = data.iloc[:,-1]

# 标准化自变量
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 划分训练集和测试集
X_train_scaled, X_test_scaled, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

# 构建神经网络模型
nn_model = MLPRegressor(hidden_layer_sizes=(64,64), max_iter=500, random_state=42)

# 训练模型
nn_model.fit(X_train_scaled, y_train)

# 测试集预测
y_pred_nn = nn_model.predict(X_test_scaled)

# 评估模型性能
mse_nn = mean_squared_error(y_test, y_pred_nn)
r2_nn = r2_score(y_test, y_pred_nn)
print(f"均方误差(MSE): {mse_nn}")
print(f"决定系数 (R²): {r2_nn}")

# 预测
new_data = [[0.84, 0.7, 0.5, 0.7]]
new_data_scaled = scaler.transform(new_data)
predict_value = nn_model.predict(new_data_scaled)
print(f"预测值：{predict_value}")

import matplotlib.pyplot as plt
import numpy as np

# 替换为实际的测试集和预测结果
y_test = np.random.uniform(3.8, 4.2, 100)  # 示例数据，替换为实际 y_test
y_pred_nn = y_test + np.random.normal(0, 0.05, 100)  # 示例数据，替换为实际预测结果

# 1. 实际值 vs. 预测值
plt.figure(figsize=(8, 6))
plt.scatter(y_test, y_pred_nn, alpha=0.7, color='blue', label='y_pred')
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', label='y_test', linewidth=2)
plt.title('y_test vs. y_pred', fontsize=16)
plt.xlabel('y_test', fontsize=14)
plt.ylabel('y_pred', fontsize=14)
plt.legend(loc='upper left', fontsize=12)
# plt.grid()
plt.show()

# 3. 神经网络训练损失曲线
# 如果使用 MLPRegressor，loss_curve_ 可用，否则模拟数据如下
loss_curve = np.logspace(0, -1, 50)  # 示例损失曲线，替换为 nn_model.loss_curve_
plt.figure(figsize=(8, 6))
plt.plot(loss_curve, 'b-', label='损失值')
plt.title('神经网络训练损失曲线', fontsize=16)
plt.xlabel('迭代次数', fontsize=14)
plt.ylabel('损失值', fontsize=14)
plt.legend(loc='upper right', fontsize=12)
plt.grid()
plt.show()


