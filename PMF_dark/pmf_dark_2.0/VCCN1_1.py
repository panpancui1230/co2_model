import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# voltage = np.array([-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0])
# mean_values = np.array([-14.8295, -12.99922, -11.54216, -10.57937, -9.12894, 
#                         -7.59989, -6.25137, -4.83564, -3.28527, -1.79184, -0.10539])
# voltage = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
# mean_values = np.array([-0.10539, 1.65365, 3.62434, 5.84268, 8.29092, 
#                         10.84378, 13.86322, 16.91701, 19.79518, 22.94266, 26.64795])
voltage = np.array([-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 
                    0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
mean_values = np.array([-14.8295, -12.99922, -11.54216, -10.57937, -9.12894, 
                        -7.59989, -6.25137, -4.83564, -3.28527, -1.79184, 
                        -0.10539, 1.65365, 3.62434, 5.84268, 8.29092, 
                        10.84378, 13.86322, 16.91701, 19.79518, 22.94266, 26.64795])

def exponential_fit(V, a, b, c):
    """
    M(V) = a * exp(b * V) + c
    """
    return a * np.exp(b * V) + c

voltage_normalized = voltage / 100.0

popt, pcov = curve_fit(exponential_fit, voltage_normalized, mean_values, maxfev=10000)

a, b, c = popt

fitted_values = exponential_fit(voltage_normalized, a, b, c)

plt.figure(figsize=(8, 6))
plt.scatter(voltage, mean_values, label='Data', color='blue', marker='o')
plt.plot(voltage, fitted_values, label=f'Fit: M(V) = {a:.3f} * exp({b:.3f} * (V/100)) + {c:.3f}', color='red', linewidth=2)
plt.xlabel('Voltage (mV)', fontsize=14)
plt.ylabel('Mean', fontsize=14)
plt.title('Exponential Fit of Voltage vs. Mean', fontsize=16)
plt.legend(fontsize=12)
plt.grid(False)
plt.tight_layout()
plt.show()

print("拟合参数：")
print(f"a = {a:.6f}")
print(f"b = {b:.6f}")
print(f"c = {c:.6f}")

print("\n拟合公式:")
print(f"M(V) = {a:.3f} * exp({b:.3f} * (V / 100)) + {c:.3f}")