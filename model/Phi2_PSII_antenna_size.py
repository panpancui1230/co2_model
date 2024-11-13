import pandas as pd

# Load the data
df = pd.read_csv('./result_cpp/result_100.csv')  # Replace with your actual file path

df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

# Calculate the new column
df['Phi2*PSII_antenna_size*100'] = df['Phi2'] * df['PSII_antenna_size']*df['ratio_absorb'] * 100

# Save the modified dataframe to a new CSV file
df.to_csv('./result_cpp/result_100_1.csv', index=False)

print("Calculation completed and saved to 'modified_file_with_results.csv'")