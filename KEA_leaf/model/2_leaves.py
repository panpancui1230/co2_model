import pandas as pd
import matplotlib.pyplot as plt


file_path_500 = './logs_500/500_CO2_plt.csv'
file_path_2000='./logs_2000/2000_CO2_plt.csv'
data_500 = pd.read_csv(file_path_500)
data_2000 = pd.read_csv(file_path_2000)


leaf_1_data_500 = data_500[data_500['叶片数'] == 1]
leaf_1_data_2000 = data_2000[data_2000['叶片数'] == 1]


plt.figure(figsize=(10, 6))
# plt.plot(leaf_1_data['ratio_absorb'], leaf_1_data['CO2_sum'], marker='o', linewidth=2,color='#FEB244')
# plt.plot(leaf_1_data['ratio_absorb'], leaf_1_data['CO2_sum'], marker='o', linewidth=2,color='green')

# 500
plt.plot(leaf_1_data_500['ratio_absorb'], leaf_1_data_500['CO2_sum'], 
         marker='o', linewidth=2, color='green', label='Light Intensity = 500')

# 20000
plt.plot(leaf_1_data_2000['ratio_absorb'], leaf_1_data_2000['CO2_sum'], 
         marker='s', linewidth=2, color='blue', label='Light Intensity = 2000')


plt.xlabel('ratio_absorb', fontsize=14)
# plt.ylabel('Total CO2', fontsize=14)
plt.ylabel(r'Total CO$_2$', fontsize=14) 
# plt.title('CO2_sum vs. ratio_absorb when Number of Leaves = 1', fontsize=16)
# plt.grid(alpha=0.5)
plt.legend(loc='lower right',bbox_to_anchor=(0.9, 0.2), fontsize=12, frameon=False)
# Show the plot
plt.tight_layout()
plt.show()