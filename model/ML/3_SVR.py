import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np

# 读取数据
file_path = './100_CO2.csv'
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 数据标准化
scaler_X = StandardScaler()
scaler_y = StandardScaler()

X_scaled = scaler_X.fit_transform(X)
y_scaled = scaler_y.fit_transform(y.values.reshape(-1, 1)).ravel()

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_scaled, test_size=0.2, random_state=42)

# 创建并训练 SVR 模型
svr = SVR(kernel='rbf', C=1.0, epsilon=0.1)
svr.fit(X_train, y_train)

# 预测
y_train_pred = svr.predict(X_train)
y_test_pred = svr.predict(X_test)

# 将预测结果逆缩放回原始范围
y_train_pred_orig = scaler_y.inverse_transform(y_train_pred.reshape(-1, 1)).ravel()
y_test_pred_orig = scaler_y.inverse_transform(y_test_pred.reshape(-1, 1)).ravel()
y_train_orig = scaler_y.inverse_transform(y_train.reshape(-1, 1)).ravel()
y_test_orig = scaler_y.inverse_transform(y_test.reshape(-1, 1)).ravel()

# 评估模型
train_mse = mean_squared_error(y_train_orig, y_train_pred_orig)
test_mse = mean_squared_error(y_test_orig, y_test_pred_orig)
train_r2 = r2_score(y_train_orig, y_train_pred_orig)
test_r2 = r2_score(y_test_orig, y_test_pred_orig)

# 输出结果
print(f"训练集均方误差 (MSE): {train_mse:.4f}")
print(f"测试集均方误差 (MSE): {test_mse:.4f}")
print(f"训练集 R²: {train_r2:.4f}")
print(f"测试集 R²: {test_r2:.4f}")

