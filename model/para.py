import pandas as pd
import itertools

# Define each variable's values
ratio_absorb_values = [0.1, 0.1889, 0.2778, 0.3667, 0.4556, 0.5444, 0.6333, 0.7222, 0.8111, 0.9]
PSII_content_per_leaf_values = [0.07, 0.1, 0.25, 0.4, 0.55, 0.7, 1, 3, 5, 7]
PSII_antenna_size_values = [0.1, 0.1889, 0.2778, 0.3667, 0.4556, 0.5444, 0.6333, 0.7222, 0.8111, 0.9]
P700_red_initial_values = [0.1, 0.4, 0.7, 1.0, 1.5, 2, 4, 6, 8, 10]

# Generate combinations
combinations = list(itertools.product(ratio_absorb_values, PSII_content_per_leaf_values, PSII_antenna_size_values, P700_red_initial_values))

# Create DataFrame
df_combinations = pd.DataFrame(combinations, columns=['ratio_absorb', 'PSII_content_per_leaf', 'PSII_antenna_size', 'P700_red_initial'])

# Save to CSV
df_combinations.to_csv("./data/variable_combinations_10000.csv", index=False)