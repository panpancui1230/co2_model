import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


data_path = "/Users/panpan/Desktop/co2_model/model/ML/500_CO2.csv"  # 替换为你的文件路径
df = pd.read_csv(data_path)


columns_to_plot = ['ratio_absorb', 'PSII_content_per_leaf', 
                   'PSII_antenna_size', 'P700_red_initial', 'CO2']
df = df[columns_to_plot]


pairplot = sns.pairplot(
    df, 
    diag_kind="hist", 
    corner=False,     
    plot_kws={'alpha': 0.5}  
)

# plt.figure(figsize=(8, 6))
plt.subplots_adjust(top=0.95)


# output_path = "pairplot_output.png" 
# pairplot.savefig(output_path)


plt.show()