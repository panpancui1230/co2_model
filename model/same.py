import pandas as pd

# Load the uploaded file
file_path_1 = './logs_cpp/combined_100_simulated_cpp_try.csv'
df1 = pd.read_csv(file_path_1)

# Load the second file (replace with actual path)
file_path_2 = './logs_QA/combined_100_simulated.csv'  # Replace with your second file's path
df2 = pd.read_csv(file_path_2)

# Ensure both files have the same columns for comparison
if set(df1.columns) != set(df2.columns):
    print("The files have different columns and cannot be compared directly.")
else:
    # Compare the data values in each column
    comparison_results = {}
    for column in df1.columns:
        comparison = df1[column].equals(df2[column])
        comparison_results[column] = "Match" if comparison else "Difference Found"
    
    # Print the comparison results
    for column, result in comparison_results.items():
        print(f"Column '{column}': {result}")
    
    # Optionally, display rows where differences were found
    differences = df1[df1 != df2]
    differences = differences.dropna(how='all')  # Remove rows with no differences
    if not differences.empty:
        print("Differences found in the following rows:")
        print(differences)
    else:
        print("All values match across both files.")