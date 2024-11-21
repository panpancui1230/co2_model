import pandas as pd

# 定义补充计算的函数
def calculate_lights_until_below_one(base_light, ratio):
    """
    计算每一层的光强,直到光强小于1为止。
    """
    results = []
    current_light = base_light * (1 - ratio)  # 第一层光强
    layer = 1
    while current_light >= 1:  # 当 light 值大于等于 1 时继续计算
        results.append({"ratio_absorb": ratio, "light": current_light, "层": layer})
        current_light *= (1 - ratio)  # 按吸收比例递减
        layer += 1
    return results

# 加载用户上传的原始数据
input_file = './light.csv'
df = pd.read_csv(input_file)

# 定义初始条件
initial_light = 2000  # 初始光强
ratios_to_calculate = [0.4, 0.3, 0.2, 0.1]  # 待补充的 ratio_absorb 值

# 计算新的数据
new_data = []
for ratio in ratios_to_calculate:
    new_data.extend(calculate_lights_until_below_one(initial_light, ratio))

# 将新数据转换为 DataFrame
new_df = pd.DataFrame(new_data)

# 合并原始数据和新数据
combined_df = pd.concat([df, new_df], ignore_index=True)

# 按照 ratio_absorb 和 层 排序
combined_df = combined_df.sort_values(by=["ratio_absorb", "层"]).reset_index(drop=True)

# 保存更新后的文件
output_file = './light_updated_v1.csv'
combined_df.to_csv(output_file, index=False)

print(f"更新后的文件已保存至: {output_file}")