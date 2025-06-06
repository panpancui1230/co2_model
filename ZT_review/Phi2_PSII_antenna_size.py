import pandas as pd

# Load the data
df = pd.read_csv('./WT_500uE_simulated.csv')  # Replace with your actual file path

# df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

# Calculate the new column
# df['Phi2*PSII_antenna_size*100'] = df['Phi2'] * df['PSII_antenna_size'] * 100
# df['Phi2*PSII_antenna_size*ratio_absorb*500'] = df['Phi2'] * 0.524 * 0.93 * 500
# df['NADPH']=df['Phi2*PSII_antenna_size*ratio_absorb*500']/2
# df['2NADPH']=df['NADPH']/2


df['3_v_ATP_made']=df['v_ATP_made']/3
# df = df.drop(columns=['2_v_ATP_made'])
# print("Column '2_v_ATP_made' has been removed.")
# # Save the modified dataframe to a new CSV file
df.to_csv('./WT_500uE_simulated.csv', index=False)

print("Calculation completed and saved to 'modified_file_with_results.csv'")