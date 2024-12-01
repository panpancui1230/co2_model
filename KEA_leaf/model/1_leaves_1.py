import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
# file_path = './logs_2000/2000_CO2_plt.csv'
file_path = './logs_500/500_CO2_plt.csv'
data = pd.read_csv(file_path)

# 筛选数据
filtered_data = data[data['叶片数'] <= 35]

# 获取唯一的 ratio_absorb 值
unique_ratios = filtered_data['ratio_absorb'].unique()

# 定义新的调色板
rocket_palette = sns.blend_palette(
    [
        # '#fffa61',  # 黄色
        # '#eb8e3d',
        # '#edad9d',
        # '#eba26e',
        # '#f0cf13',
        # '#d46833',
        '#f0f06e',
        # '#FFF95D',
        # '#979c70',
        '#c5d420',  # 浅绿色
        '#8ab81f',  # 明绿色
        '#1fb861',  # 青绿色
        # '#0d631a',
        '#1fb8b0',  # 蓝色
        '#1f8ab8',  # 紫色
        '#1f6eb8',  # 粉色
        # '#130969'   # 红色

    # '#ffda05',  # 鲜红
    # '#ffff05',  # 橙黄
    # '#daff05',  # 浅蓝绿
    # '#8aff05',  # 蓝
    # '#05ff54',  # 深蓝
    # '#05ffd1',  # 绿松石
    # '#05d1ff',  # 烧橙
    # '#0586ff',  # 墨蓝
    # '#0509ff'   # 白绿
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
plt.rcParams['font.family'] = 'Arial'

# 创建绘图
plt.figure(figsize=(10, 6))

# 绘制每条曲线
for ratio in unique_ratios:
    # 筛选当前 ratio_absorb 的数据
    ratio_data = filtered_data[filtered_data['ratio_absorb'] == ratio]
    
    # 绘制曲线
    plt.plot(ratio_data['叶片数'], ratio_data['CO2_sum'], 
             label=f'ratio_absorb={ratio:.1f}', 
             marker='o', markersize=5, linewidth=2.5, color=colors[ratio])

# 设置坐标轴范围和标签
plt.xlim(0, 37)
plt.xlabel('Number of overlapping leaves', fontsize=14)
plt.ylabel(r'$A_{total}$ (µmol $m^{-2} s^{-1}$)', fontsize=14)

# 添加图例
plt.legend(loc='upper left', bbox_to_anchor=(0.70, 0.55), frameon=False, fontsize=12, title_fontsize=14)

# 调整布局
plt.tight_layout()

# 显示图像
plt.show()