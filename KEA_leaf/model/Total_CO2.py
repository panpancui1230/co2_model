import pandas as pd

# 读取CSV文件
input_file = "./logs_500/2000_CO2.csv"  # 替换为实际路径
data = pd.read_csv(input_file)

# 确保文件中包含 'ratio_absorb' 和 'CO2' 列
if 'ratio_absorb' not in data.columns or 'CO2' not in data.columns:
    raise ValueError("输入文件中必须包含 'ratio_absorb' 和 'CO2' 列")

# 定义需要计算的 ratio_absorb 值
ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# 创建一个空字典存储每个 ratio_absorb 的 CO2 总和
result_dict = {"ratio_absorb": [], "CO2_sum": []}

# 遍历每个 ratio_absorb，计算对应的 CO2 总和
for ratio in ratios:
    co2_sum = data[data["ratio_absorb"] == ratio]["CO2"].sum()
    result_dict["ratio_absorb"].append(ratio)
    result_dict["CO2_sum"].append(co2_sum)

# 转换为DataFrame
result_df = pd.DataFrame(result_dict)

# 输出结果到CSV文件
output_file = "./logs_500/CO2_sum_by_ratio_absorb.csv"
result_df.to_csv(output_file, index=False)

print(f"计算完成，结果已保存到文件: {output_file}")