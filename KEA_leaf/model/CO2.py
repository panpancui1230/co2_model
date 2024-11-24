import pandas as pd

# Load the data
file_path = './logs_500/2000_ATP_made.csv'  # Replace with your actual file path
df = pd.read_csv(file_path)

# Calculate CO2 as the smaller value between 2NADPH and 2_v_ATP_made
df['CO2'] = df[['2NADPH', '3_v_ATP_made']].min(axis=1)


# Save the modified dataframe to a new file
output_path = './logs_500/2000_CO2.csv'  # Replace with your desired output file path
df.to_csv(output_path, index=False)

print(f"CO2 calculated and saved to '{output_path}'.")