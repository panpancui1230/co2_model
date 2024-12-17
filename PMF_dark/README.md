# KEA3项目结构

```bash
├── README.md
├── model（代码实现）
│   ├── KEA3_Q_ref.py
│   ├── KEA3_Q.py
│   ├── calc.py
│   ├── data（运算需要的数据）
│   │   ├── constants.csv
│   │   ├── initial_states.csv
│   ├── logs（输出的运算结果数据文件）
│   │   ├── WT_1000uE_simulated.csv
│   │   ├── WT_100uE_simulated.csv
│   │   ├── WT_250uE_simulated.csv
│   │   ├── WT_500uE_simulated.csv
│   │   ├── WT_50uE_simulated.csv
│   │   ├── delta_NPQ_LEF1000_uE_simulated.csv
│   │   ├── delta_NPQ_LEF100_uE_simulated.csv
│   │   ├── delta_NPQ_LEF250_uE_simulated.csv
│   │   ├── delta_NPQ_LEF500_uE_simulated.csv
│   │   ├── delta_NPQ_LEF50_uE_simulated.csv
│   │   ├── kea3_1000uE_simulated.csv
│   │   ├── kea3_100uE_simulated.csv
│   │   ├── kea3_250uE_simulated.csv
│   │   ├── kea3_500uE_simulated.csv
│   │   └── kea3_50uE_simulated.csv
│   ├── painter.py
│   ├── sun_sim.py
│   └── utils.py
├── paper（论文相关资料）
│   ├── Impact of ion fluxes across thylakoid membranes on photosynthetic electron transport and photoprotection..pdf
│   └── Li_2021_Nplants_SI.pdf
├── parts（几个部分代码的拆解，仅供阅读，不要运行）
│   ├── 1_Code related to naming states and constants.py
│   ├── 2_Code related to ODE simulations.py
│   ├── 3_Code related to plotting out results.py
│   ├── 4_ Initial concentrations and parameters.py
│   ├── 5_Classes to hold constants and states.py
│   ├── 6_ Display notes, states and constants.py
│   └── 7_Startup code to get things set up.py
└── tools（查看原始总代码结构的工具）
    ├── arc.md
    ├── arc.pdf
    └── arc_ana.py

```

## 1. 查看代码架构

通过[工具文档](./tools/arc.md)，可以看到当前修改代码的函数定义、类及其方法的定义，有助于对代码进行理解。

## 2. 配置环境

```bash
pip install numpy, matplotlib, scipy, pandas, copy
```

## 3. 运行代码

运行优化版：

```bash
cd ./model
python KEA3_Q.py
```

运行原始版：

```bash
cd ./model
python KEA3_Q_ref.py
```

**注意**：需要比对两个代码输出的图像，确保结果计算的正确。

## 4. 当前开发进度

- [x] 完成原版代码的冗余剔除
剔除不需要的函数、变量、类
- [x] 完成修改版代码的IO处理
实现`./data/*.CSV`文件中定义数据常量，并保存输出结果到`./log/*.CSV`文件中
- [x] 完成修改版代码的函数到类方法的映射
其中：

- 在`utils.py`中指定了数据输入输出方式及内容、初始状态和常量的设定

- 在`painter.py`中实现图像的绘制

- 在`calc.py`中实现反应需要的计算单元

- 在`sun_sim.py`中实现模拟光照的过程
