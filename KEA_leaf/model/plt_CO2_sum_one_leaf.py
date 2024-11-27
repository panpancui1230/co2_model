import pandas as pd
import matplotlib.pyplot as plt

# Load the uploaded CSV file
file_path = './logs_500/500_CO2_plt.csv'
data = pd.read_csv(file_path)

# Filter data for '叶片数' == 1
leaf_1_data = data[data['叶片数'] == 1]

# Create a line plot for CO2_sum vs. ratio_absorb when 叶片数 == 1
plt.figure(figsize=(10, 6))
# plt.plot(leaf_1_data['ratio_absorb'], leaf_1_data['CO2_sum'], marker='o', linewidth=2,color='#FEB244')
plt.plot(leaf_1_data['ratio_absorb'], leaf_1_data['CO2_sum'], marker='o', linewidth=2,color='green')

# Add labels, title, and grid
plt.xlabel('ratio_absorb', fontsize=14)
plt.ylabel('Total CO2', fontsize=14)
# plt.title('CO2_sum vs. ratio_absorb when Number of Leaves = 1', fontsize=16)
# plt.grid(alpha=0.5)

# Show the plot
plt.tight_layout()
plt.show()