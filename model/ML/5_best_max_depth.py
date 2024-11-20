import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

# 读取数据
file_path = './100_CO2.csv'
df = pd.read_csv(file_path)

# 分割特征 (X) 和目标变量 (y)
X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]

# 数据集划分 (80% 训练集, 20% 测试集)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 使用 GridSearchCV 寻找最佳 max_depth
param_grid = {'max_depth': [5, 10, 15, 20, None]}
grid_search = GridSearchCV(RandomForestRegressor(random_state=42), param_grid, cv=5, scoring='r2')
grid_search.fit(X_train, y_train)

# 输出最佳参数和对应的分数
best_max_depth = grid_search.best_params_['max_depth']
best_r2_score = grid_search.best_score_
print(f"Best max_depth: {best_max_depth}")
print(f"Best R² Score from GridSearchCV: {best_r2_score:.4f}")

# 使用最佳参数训练模型
rf = RandomForestRegressor(n_estimators=100, max_depth=best_max_depth, random_state=42)
rf.fit(X_train, y_train)

# 模型预测
y_train_pred = rf.predict(X_train)  # 训练集预测
y_test_pred = rf.predict(X_test)   # 测试集预测

# 训练集评估
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