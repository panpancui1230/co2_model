import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error


file_path = './500_CO2.csv'
df = pd.read_csv(file_path)

plt.rcParams['font.sans-serif'] ='Arial'
plt.rcParams["font.size"] = 18

X = df[["ratio_absorb", "PSII_content_per_leaf", "PSII_antenna_size", "P700_red_initial"]]
y = df["CO2"]


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


param_grid = {
    "n_estimators": [420,450,480],
    "max_depth": [19,20, 21],
    "min_samples_split": [2, 3],
    "min_samples_leaf": [1, 2],
    "max_features": ["sqrt", "log2"]
}


grid_search = GridSearchCV(
    estimator=RandomForestRegressor(random_state=42),
    param_grid=param_grid,
    cv=5,
    scoring="neg_mean_squared_error",
    n_jobs=-1
)


grid_search.fit(X_train, y_train)

# 最佳参数
best_params = grid_search.best_params_
print("Best Parameters:", best_params)

# 使用最佳参数训练模型
optimized_rf = grid_search.best_estimator_

# 预测
y_train_pred = optimized_rf.predict(X_train)
y_test_pred = optimized_rf.predict(X_test)

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
print(f"Train R-squared (R2): {train_r2:.4f}")

print("\nTesting Data Evaluation:")
print(f"Test Mean Squared Error (MSE): {test_mse:.4f}")
print(f"Test Root Mean Squared Error (RMSE): {test_rmse:.4f}")
print(f"Test Mean Absolute Error (MAE): {test_mae:.4f}")
print(f"Test R-squared (R2): {test_r2:.4f}")

# 残差图
residuals = y_test - y_test_pred
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(8, 6))
sns.scatterplot(x=y_test_pred, y=residuals, color='green', alpha=0.7)
plt.axhline(0, color='red', linestyle='--')
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.title("Residual Plot")
plt.show()

# 实际值 vs 预测值
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(8, 6))
sns.scatterplot(x=y_test, y=y_test_pred, color='blue', alpha=0.7)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', linewidth=2)
plt.xlabel("Actual A Values (µmol $m^{-2} s^{-1}$)")
plt.ylabel("Predicted A Values (µmol $m^{-2} s^{-1}$)")
plt.title("Actual vs Predicted Values")
plt.show()

# MSE
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


