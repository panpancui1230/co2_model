import pandas as pd

# Load the data
df = pd.read_csv('./data/result_7.0.csv')  # Replace with your actual file path

# Calculate the new column
df['Result'] = df['Phi2'] * df['PSII_antenna_size'] * 100

# Save the modified dataframe to a new CSV file
df.to_csv('./data/result_7.0.csv', index=False)

print("Calculation completed and saved to 'modified_file_with_results.csv'")