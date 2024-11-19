import pandas as pd

file_path = './data/variable_combinations_idx.csv'
df = pd.read_csv(file_path)

# 添加 idx 列，从 0 开始
# df["idx"] = range(len(df))

#删除
# df_filtered = df[~df["idx"].isin([1093,2017,3007,3043,4017,5039,6049,6058,7013,
#                                   8049,8248,9022,9033,9059,9137,9237])]

df_filtered = df[~df["idx"].isin([69,3268,3573,4335,5026,5218,5341,5348,5482,6012,6154,
                                  6762,7233,7418,7555,7639,8008,8237,8757,9445])]
# 如果需要保存为 CSV
df_filtered.to_csv('./data/variable_combinations_idx_500.csv', index=False)