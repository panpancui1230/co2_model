
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = './logs_2/2000_CO2_plt.csv'
data = pd.read_csv(file_path)

# Extract unique values for `ratio_absorb`
unique_ratios = data['ratio_absorb'].unique()

# Define a color map, where we swap the colors for ratio=0.1 and ratio=0.3
colors = {}
default_colors = plt.cm.tab10(range(len(unique_ratios)))
for idx, ratio in enumerate(unique_ratios):
    if ratio == 0.1:
        colors[ratio] = 'tab:green'  # Swap to another specific color
    elif ratio == 0.3:
        colors[ratio] = 'tab:blue'   # Swap to another specific color
    else:
        colors[ratio] = default_colors[idx % len(default_colors)]

# Create a line plot for CO2_sum against 叶片数 for each unique ratio_absorb value
plt.figure(figsize=(12, 6))

for ratio in unique_ratios:
    # Filter data for the current ratio_absorb
    filtered_data = data[data['ratio_absorb'] == ratio]
    
    # Plot the line and connect the points
    plt.plot(filtered_data['叶片数'], filtered_data['CO2_sum'], 
             label=f'ratio_absorb={ratio:.2f}', 
             marker='o',markersize=4,linewidth=2, color=colors[ratio])

# Add labels, title, legend, and grid
plt.xlabel('Number of leaves', fontsize=14)
plt.ylabel('Total CO2', fontsize=14)
# plt.title('CO2_sum vs. Number of leaves for Different ratio_absorb Values', fontsize=16)
plt.legend(title='Ratio Absorb', loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

# Show the plot
plt.show()