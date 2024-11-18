import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error,r2_score

# 加载数据
file_path = '/Users/panpan/Desktop/co2_model/model/ML/result_100.csv'
data = pd.read_csv(file_path)

# 自变量和因变量
X = data.iloc[:,:-1]
y = data.iloc[:,-1]

# 标准化自变量
scaler = RobustScaler()
X_scaled = scaler.fit_transform(X)

# 划分训练集和测试集
poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform
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

