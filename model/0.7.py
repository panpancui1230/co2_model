import pandas as pd

# Load the CSV file
file_path = './data/initial_states.csv'
df = pd.read_csv(file_path)

# Update the 'P700_red_initial' column to 0.7 for all rows
df['P700_red_initial'] = 0.7

# Save the modified DataFrame to a new CSV file
updated_file_path = './data/initial_states_updated_P700_red.csv'
df.to_csv(updated_file_path, index=False)

updated_file_path