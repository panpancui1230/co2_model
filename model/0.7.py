# import pandas as pd

# # Load the CSV file
# file_path = './data/result_2.0.csv'
# df = pd.read_csv(file_path)

# # Update the 'P700_red_initial' column to 0.7 for all rows
# df['P700_red_initial'] = 0.7

# # Save the modified DataFrame to a new CSV file
# updated_file_path = './data/result_2.0.csv'
# df.to_csv(updated_file_path, index=False)

# updated_file_path



# import pandas as pd

# # Load the data
# df = pd.read_csv('./data/initial_states_6.0.csv')  # Replace with the actual file path

# # Update values in 'P700_red_initial' column
# df['P700_red_initial'] = df['P700_red_initial'].replace({0.1:0.7, 0.5:1, 1:2})

# # Save the modified dataframe to a new CSV file
# df.to_csv('./data/initial_states_7.0.csv', index=False)

# print("Values updated and saved to 'updated_file.csv'")



import pandas as pd

# 读取数据
df = pd.read_csv('./data/initial_states_6.0.csv')  # 确保路径正确

# 修改最后 100 行的 P700_red_initial 列的值为 0.1
df.loc[-100:, 'P700_red_initial'] = 0.7

# 将修改后的数据保存到新的 CSV 文件
df.to_csv('./data/initial_states_7.0.csv', index=False)

print("最后 100 个 P700_red_initial 值已更新为 0.1,并保存到 'initial_states_7.0.csv'")