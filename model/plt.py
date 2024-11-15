import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = './logs_QA/result_100_CO2.csv'  # Replace with your actual file path
df = pd.read_csv(file_path)

# Check if required columns are present
if 'idx' in df.columns and 'CO2' in df.columns:
    # Plot CO2 vs idx
    plt.figure(figsize=(8, 5))
    plt.plot(df['idx'], df['CO2'], marker='o', linestyle='-', label='CO2')
    plt.title("CO2 vs idx", fontsize=14)
    plt.xlabel("idx", fontsize=12)
    plt.ylabel("CO2", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print("The file does not contain the required columns 'idx' and 'CO2'.")