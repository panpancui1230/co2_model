import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = './logs_2/2000_CO2_plt.csv'
data = pd.read_csv(file_path)

# Extract unique values for `ratio_absorb`
unique_ratios = data['ratio_absorb'].unique()

# Plot CO2_sum against 叶片数 for each unique ratio_absorb value
plt.figure(figsize=(12, 6))

for ratio in unique_ratios:
    # Filter data for the current ratio_absorb
    filtered_data = data[data['ratio_absorb'] == ratio]
    
    # Plot the data
    plt.plot(filtered_data['叶片数'], filtered_data['CO2_sum'], label=f'ratio_absorb={ratio:.2f}',linewidth=3)

# Add labels, title, legend, and grid
plt.xlabel('Number of leaves', fontsize=14)
plt.ylabel('CO2_sum', fontsize=14)
plt.title('CO2_sum vs. Number of leaves for Different ratio_absorb Values', fontsize=16)
plt.legend(title='Ratio Absorb', loc='upper left', bbox_to_anchor=(1, 1))
# plt.grid(alpha=0.5)
plt.tight_layout()

# Show the plot
plt.show()