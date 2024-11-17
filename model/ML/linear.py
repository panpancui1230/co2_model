import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# 加载数据
file_path = '/Users/panpan/Desktop/ML/result_100.csv'

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

# 输出性能评估指标
print(f"均方误差 (MSE): {mse}")
print(f"决定系数 (R²): {r2}")

# 查看模型的系数和截距
coefficients = lr_model.coef_  # 每个自变量的系数
intercept = lr_model.intercept_  # 截距

print("模型系数:", coefficients)
print("模型截距:", intercept)

# 预测新数据
new_data = [[0.84, 0.7, 0.5, 0.7]] 
predicted_value = lr_model.predict(new_data)
print(f"新数据的预测值: {predicted_value[0]}")