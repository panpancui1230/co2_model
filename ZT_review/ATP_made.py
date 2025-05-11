import pandas as pd

# Load the CSV file
file_path = './WT_500uE_simulated.csv'
df = pd.read_csv(file_path)

# Define the function to calculate v_ATP_made
def Vproton_pmf_actvt(pmf):
    # Fraction of activity based on pmf
    v_proton_active = 1 - (1 / (10 ** ((pmf - 0.132) * 1.5 / 0.06) + 1))  # Reduced ATP synthase
    v_active = v_proton_active * 4.667 * 200
    v_proton_ATP = v_active
    v_ATP_made = v_proton_ATP / 4.667
    return v_ATP_made

# Ensure 'pmf' column exists
if 'pmf' not in df.columns:
    raise ValueError("The column 'pmf' does not exist in the input file.")

# Calculate v_ATP_made and store it in a new column
df['v_ATP_made'] = df['pmf'].apply(Vproton_pmf_actvt)

# Save the modified dataframe to a new CSV file
output_path = './WT_500uE_simulated.csv'
df.to_csv(output_path, index=False)

print(f"v_ATP_made calculated and saved to '{output_path}'.")