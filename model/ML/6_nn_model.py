import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error,r2_score, mean_absolute_error
import matplotlib.pyplot as plt
import numpy as np

file_path = './500_CO2.csv'
data = pd.read_csv(file_path)
plt.rcParams['font.sans-serif'] ='Arial'
plt.rcParams["font.size"] = 18

X = data.iloc[:,:-1]
y = data.iloc[:,-1]


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 构建神经网络模型
# nn_model = MLPRegressor(hidden_layer_sizes=(64,64), max_iter=500, random_state=42)
nn_model = MLPRegressor(hidden_layer_sizes=(256,128,64), max_iter=50,random_state=42)
# nn_model = MLPRegressor(
#     hidden_layer_sizes=(256, 128, 64),
#     max_iter=50,
#     alpha=0.01,  # L2 正则化
#     learning_rate="adaptive",
#     learning_rate_init=0.01,
#     early_stopping=True,
#     validation_fraction=0.2,
#     random_state=42
# )


nn_model.fit(X_train_scaled, y_train)

# 测试集预测
y_train_pred = nn_model.predict(X_train_scaled)

y_test_pred = nn_model.predict(X_test_scaled)

# 训练集评估
train_mse = mean_squared_error(y_train, y_train_pred)
train_rmse = train_mse ** 0.5
train_mae = mean_absolute_error(y_train, y_train_pred)
train_r2 = r2_score(y_train, y_train_pred)

# 测试集评估
test_mse = mean_squared_error(y_test, y_test_pred)
test_rmse = test_mse ** 0.5
test_mae = mean_absolute_error(y_test, y_test_pred)
test_r2 = r2_score(y_test, y_test_pred)


print("\nTraining Data Evaluation:")
print(f"Train Mean Squared Error (MSE): {train_mse:.4f}")
print(f"Train Root Mean Squared Error (RMSE): {train_rmse:.4f}")
print(f"Train Mean Absolute Error (MAE): {train_mae:.4f}")
print(f"Train R-squared (R²): {train_r2:.4f}")

print("\nTesting Data Evaluation:")
print(f"Test Mean Squared Error (MSE): {test_mse:.4f}")
print(f"Test Root Mean Squared Error (RMSE): {test_rmse:.4f}")
print(f"Test Mean Absolute Error (MAE): {test_mae:.4f}")
print(f"Test R-squared (R²): {test_r2:.4f}")




#残差图 (Residual Plot)

residuals = y_test - y_test_pred


residuals_df = pd.DataFrame({
    "Predicted Values": y_test_pred,
    "Residuals": residuals
})

plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(8, 6))
sns.scatterplot(data=residuals_df, x="Predicted Values", y="Residuals", color='green', alpha=0.7, label="Residuals")
plt.axhline(y=0, color='r', linestyle='--', label="Residual = 0") 

plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.legend(loc='upper left')
plt.grid(False)
plt.show()

#实际值 vs 预测值散点图
results_df = pd.DataFrame({
    "Actual": y_test,
    "Predicted": y_test_pred
})

plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(8, 6))
sns.scatterplot(data=results_df, x="Actual", y="Predicted", color='blue', alpha=0.7, label="Predicted vs Actual")
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', label="Ideal Prediction")  # Ideal prediction line


plt.xlabel("Actual A Values (µmol $m^{-2} s^{-1}$)")
plt.ylabel("Predicted A Values (µmol $m^{-2} s^{-1}$)")
plt.title("Predicted vs Actual Values")
plt.legend(loc='upper left')
plt.grid(False)
plt.show()

# 3. 神经网络训练损失曲线
plt.rcParams['font.family'] = 'Arial'
loss_curve = np.logspace(0, -1, 50)  
plt.figure(figsize=(8, 6))
plt.plot(loss_curve, 'b-', label='loss')
plt.title('loss curve', fontsize=16)
plt.xlabel('epoch', fontsize=14)
plt.ylabel('loss', fontsize=14)
plt.legend(loc='upper right', fontsize=12)
plt.grid(False)
plt.show()

#metrics
plt.rcParams['font.family'] = 'Arial'
metrics = ["MSE", "RMSE", "MAE", "R^2"]
train_values = [train_mse, train_rmse, train_mae, train_r2]
test_values = [test_mse, test_rmse, test_mae, test_r2]


plt.figure(figsize=(10, 6))
plt.plot(metrics, train_values, marker='o', label="Training Data", linestyle='-', linewidth=2, markersize=8, color='blue')
plt.plot(metrics, test_values, marker='o', label="Testing Data", linestyle='-', linewidth=2, markersize=8, color='green')


for i, value in enumerate(train_values):
    plt.text(i, value + 0.01, f"{value:.4f}", ha="center", va="bottom", fontsize=10, color='blue')
for i, value in enumerate(test_values):
    plt.text(i, value + 0.01, f"{value:.4f}", ha="center", va="bottom", fontsize=10, color='green')


plt.title("Model Evaluation Metrics: Training vs Testing", fontsize=16)
plt.ylabel("Metric Value", fontsize=14)
plt.xlabel("Metrics", fontsize=14)


plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()


