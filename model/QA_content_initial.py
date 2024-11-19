import pandas as pd

# Load the data
df = pd.read_csv('./data/variable_combinations.csv')

#QA保留3位小数
df['QA_content_initial'] = (df['PSII_content_per_leaf']/0.7).round(5)

# Save
df.to_csv('./data/variable_combinations.csv', index=False)

print("Values updated and saved to 'updated_file.csv'")


