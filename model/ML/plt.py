import matplotlib.pyplot as plt
import numpy as np

# 替换为实际的测试集和预测结果
y_test = np.random.uniform(3.8, 4.2, 100)  # 示例数据，替换为实际 y_test
y_pred_nn = y_test + np.random.normal(0, 0.05, 100)  # 示例数据，替换为实际预测结果

# 1. 实际值 vs. 预测值
plt.figure(figsize=(8, 6))
plt.scatter(y_test, y_pred_nn, alpha=0.7, color='blue', label='预测点')
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', label='理想预测', linewidth=2)
plt.title('实际值 vs. 预测值', fontsize=16)
plt.xlabel('实际值 (y_test)', fontsize=14)
plt.ylabel('预测值 (y_pred)', fontsize=14)
plt.legend(loc='upper left', fontsize=12)
plt.grid()
plt.show()
