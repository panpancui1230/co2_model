import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
# file_path = './logs_2000/2000_CO2_plt.csv'
file_path = './logs_500/500_CO2_plt.csv'
data = pd.read_csv(file_path)
plt.rcParams['font.sans-serif'] ='Arial'
plt.rcParams["font.size"] = 18
# 筛选数据
# filtered_data = data[data['叶片数'] <= 35]
filtered_data = data[data['叶片数'] <= 10]

# 获取唯一的 ratio_absorb 值
unique_ratios = filtered_data['ratio_absorb'].unique()


rocket_palette = sns.blend_palette(
    [
        
        '#2E5FAD',  # 更深的蓝色
        '#3FA653',  # 更深的绿色
        '#E4B627',  # 更深的黄色
        '#D14B26'   # 更深的橙色
    ],
    n_colors=len(unique_ratios)
)

# 根据新的调色板生成颜色字典
colors = {ratio: rocket_palette[idx % len(rocket_palette)] for idx, ratio in enumerate(unique_ratios)}
# 使用 YlGn 调色板生成调色板，从黄到绿渐变
# ylgn_palette=sns.color_palette("rocket_r", as_cmap=True)
# ylgn_palette = sns.color_palette("YlGn", n_colors=len(unique_ratios))
# ylgn_palette = sns.dark_palette("green", n_colors=len(unique_ratios), reverse=True)
# ylgn_palette=sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, as_cmap=True)
# colors = {ratio: ylgn_palette[idx % len(ylgn_palette)] for idx, ratio in enumerate(unique_ratios)}

# rocket_palette = sns.color_palette("plasma", n_colors=len(unique_ratios))

# rocket_palette = sns.blend_palette(["#f16c23",'#eba95e',"#1b7c3d","#2b6a99"], n_colors=len(unique_ratios))
# rocket_palette = sns.blend_palette(["#fe2f2f",'#fe7777',"#9797f8","#2a2afc"], n_colors=len(unique_ratios))

# rocket_palette = sns.blend_palette(["#D4E157","#DCE775","#8BC34A","#7CB342","#388E3C","#1B5E20"], n_colors=len(unique_ratios))
# rocket_palette = sns.blend_palette(["#FFEB3B","#FFF176","#FFF59D","#DCE775","#CDDC39","#388E3C","#1B5E20"], n_colors=len(unique_ratios)) 

# rocket_palette = sns.blend_palette(['#FEDACC', '#FD9778', '#D5E5F4', '#93C7E6','#C9EBC4', '#9EDB99'], n_colors=len(unique_ratios))
# rocket_palette = sns.blend_palette(['#FFF59D','#FFEE58','#C5E1A5','#8BC34A','#33691E','#4292C6','#08306B'], n_colors=len(unique_ratios))

# rocket_palette = sns.blend_palette(["#daed0c", "#0c6e10"], n_colors=len(unique_ratios))
# colors = {ratio: rocket_palette[idx % len(rocket_palette)] for idx, ratio in enumerate(unique_ratios)}
# ylgn_palette = sns.color_palette("YlGn", n_colors=len(unique_ratios))
# colors = {ratio: ylgn_palette[idx] for idx, ratio in enumerate(unique_ratios)}

# 使用 blend 调色板（从蓝绿色到橙色）
# blend_palette = sns.color_palette("blend:#7AB,#EDA", n_colors=len(unique_ratios))
# colors = {ratio: blend_palette[idx % len(blend_palette)] for idx, ratio in enumerate(unique_ratios)}


# 使用绿色渐变调色板
# green_palette = sns.color_palette("Greens", n_colors=len(unique_ratios))
# colors = {ratio: green_palette[idx % len(green_palette)] for idx, ratio in enumerate(unique_ratios)}

# ylgn_palette = sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, n_colors=len(unique_ratios))
# colors = {ratio: ylgn_palette[idx % len(ylgn_palette)] for idx, ratio in enumerate(unique_ratios)}
# 设置字体
# 设置字体
# 设置字体
plt.rcParams['font.family'] = 'Arial'

# 创建绘图
plt.figure(figsize=(10, 6))

# 定义 3 种点样式
markers = ['o', 's', 'D']  # 圆形、方形、菱形

# 绘制每条曲线
for idx, ratio in enumerate(unique_ratios):
    # 筛选当前 ratio_absorb 的数据
    ratio_data = filtered_data[filtered_data['ratio_absorb'] == ratio]
    
    # 绘制曲线
    plt.plot(ratio_data['叶片数'], ratio_data['CO2_sum'], 
            #  label=f'ratio_absorb={ratio:.1f}', 
             label=f'$Abs_{{ratio}}={ratio:.1f}$',
             marker=markers[idx % len(markers)],  # 循环使用 3 种点样式
             markersize=7.5, linewidth=2.5, color=colors[ratio])  # 保持颜色不变

# 设置坐标轴范围和标签
# plt.xlim(0, 37)
plt.xlim(0, 15)
plt.xlabel('Number of overlapping leaves', fontsize=18)
plt.ylabel(r'$A_{total}$ (µmol $m^{-2} s^{-1}$)', fontsize=18)

# 添加图例
plt.legend(loc='upper left', bbox_to_anchor=(0.72, 0.62), frameon=False, fontsize=14.5, title_fontsize=14)
# plt.legend(loc='lower right', bbox_to_anchor=(0.30, 0.62), frameon=False, fontsize=14.5, title_fontsize=14)
# 调整布局
plt.tight_layout()

# 显示图像
plt.show()