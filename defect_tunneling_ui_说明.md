# defect_tunneling_ui 使用说明

本文档对应本地可执行页面：

```text
Current_simulate\defect_tunneling_ui.html
Current_simulate\run_defect_tunneling_ui.py
```

该页面用于展示 FN 隧穿电流、缺陷俘获电荷、局域场变化、势垒形状变化和 KMC 导电通道形成之间的关系。它是 reduced-order 可视化模型，不是 Comphy 仿真，也不是完整 microscopic TAT / NMP / Poisson 自洽求解器。

## 1. 运行方式

在项目根目录运行：

```powershell
py -3.12 -S .\Current_simulate\run_defect_tunneling_ui.py
```

打开：

```text
http://127.0.0.1:8765/defect_tunneling_ui.html
```

只启动服务、不自动打开浏览器：

```powershell
py -3.12 -S .\Current_simulate\run_defect_tunneling_ui.py --no-open
```

`-S` 用于跳过用户级 Python site 初始化。本 UI 服务脚本只依赖 Python 标准库。

## 2. 页面结构

页面主区域包含四个图：

| 图 | 显示内容 | 物理含义 |
|---|---|---|
| 电流 `I(t)` 随时间变化 | 含缺陷电流、无缺陷基准、当前时间、KMC 击穿展示 | 主结果，观察电流随时间衰减和击穿跳变 |
| 氧化层厚度起伏 Monte Carlo | 多条器件 I-t、器件间 P5/P50/P95 电流包络、Weibull tBD 分布 | 展示微小 `tox` 起伏导致的器件间 FN 电流和击穿时间差异 |
| KMC 导电通道形成 | 按器件横向长度和氧化层厚度绘制的二维栅格，单元为 `1 nm × 1 nm` | 演示缺陷如何形成 gate-to-substrate 通道 |
| 局域场导致的隧穿势垒变化 | `SiO2` 导带 `Uc(x)`、氧化物场强 `Eox(x)`、缺陷点 | 观察俘获电荷如何改变局域场和势垒形状 |
| 缺陷俘获电荷与有效场 | 左轴 `Qdef(t)`，右轴 `Eeff(t)` | 解释电流衰减的中间变量 |

## 2.1 实验数据导入

页面默认加载：

```text
Current_simulate\data\previous_current.csv
```

当前默认数据按 150 °C 实验条件处理；默认加载该文件时，页面温度会设置为：

```text
T = 423.15 K
```

该文件由原始 Excel 数据：

```text
Current_simulate\1.xlsx
```

中的两列转换得到：

```text
CH11_Time_s
CH11_I_A
```

左侧 `导入实验 I-t 数据 (CSV/TXT)` 支持导入新的两列数据。支持以下格式：

- 带表头：`CH11_Time_s, CH11_I_A`
- 带通用表头：`time_s, current_a` 或类似包含 `time/current` 的列名
- 无表头：默认使用前两列作为 `time_s` 和 `I_A`

导入时会过滤非正时间、非正电流和无效数字。为避免浏览器绘图卡顿，实验数据会自动降采样到最多 2500 个点用于显示；原始 CSV 文件本身不会被修改。

导入完成后，UI 会自动执行一次非线性拟合。该拟合不再引入独立的 `K_FN` 拟合前因子，也不再通过缩放 `Area` 来做幅值匹配；`Area` 保持当前输入值不变，`Phi_B` 固定为 SiC/SiO2 导带差 `2.7 eV`，`mOx` 固定为 `0.29 m0`，`mSiC` 固定为 `0.42 m0`。自动拟合会调整以下物理/等效参数：

- `Vg`：决定氧化层平均场强 `E0 = Vg / tox`；
- 三类缺陷的 `Nsheet`：决定可俘获电荷容量；
- 各缺陷族深度 `x_k`：进入局域场屏蔽和 `tau_eff(x_k)` 计算。

各缺陷输入框中的 `tau`、`Area`、`Phi_B`、`mOx`、`mSiC`、`tox`、`epsOx`、温度、`gamma` 和 `beta` 不作为自动拟合参数；当前时间常数变化由缺陷空间深度和场强模型通过 `tau_eff` 派生，缺陷能级固定为 `Et=Ec(SiC)=0 eV` 且不参与 `tau_eff`。这样可以避免把“时间常数”当作任意数学自由度去吸收曲线误差。

目标函数为击穿前数据段上的 log 电流误差：

```text
minimize RMSE_log10 = sqrt( mean( [log10(I_exp) - log10(I_model)]^2 ) )
```

拟合使用轻量级坐标搜索，适合在浏览器中快速给出可交互的等效参数；它不是全局优化器，也不应把单次拟合结果直接解释为唯一微观缺陷参数。导入后顶部会显示两个拟合误差：

| 指标 | 含义 |
|---|---|
| `R²` | 在击穿前数据段上，用 log10 电流计算的决定系数 |
| `RMSE (ln)` | 在击穿前数据段上，`ln(I_exp) - ln(I_model)` 的均方根误差 |

误差计算会自动忽略非正电流，并在检测到电流突跳击穿后只使用击穿前的数据点，避免击穿后的导电通道电流主导 FN 衰减段拟合。

顶部指标栏只显示四个关键量：

| 指标 | 含义 |
|---|---|
| `I(t)` | 当前时间点的含缺陷 FN 电流 |
| `Eeff(t)` | 当前时间点的有效氧化层场强 |
| `KMC tBD` | KMC 贯通时间；未贯通时显示未贯通状态 |
| `KMC path` | KMC 通道状态 |

`log10 I(t)`、`I/I_no_defect`、`Qdef(t)`、`tau_eff range` 等二级信息仍在内部计算和 CSV 导出中保留，但不再占用顶部空间。

## 3. 图表刻度

所有图表采用统一的标准刻度格式：

- 时间对数轴显示为 `10^n`，不显示 `1e5`。
- 电流对数轴采用右上角倍率，例如 `×10^-5 A`，y 轴刻度显示普通数字，如 `1, 0.1, 0.01`。
- 线性轴在数值过大或过小时采用轴标题倍率，例如 `Qdef (C/m²) ×10^-3`。
- `Qdef(t)` 和 `Eeff(t)` 使用双 y 轴：左轴为 `Qdef`，右轴为 `Eeff`，避免不同量纲共用一个轴。

## 4. 基础电流模型

主电流模型为 FN-screened defect charging：

```text
I(t) = Area * A_FN * E_inj(t)^2 * exp(-B_FN / E_inj(t))
Q_def(t) = sum_k Q_k * (1 - exp(-t / tau_eff,k))
E_inj(t) = E0 - Q_def(t) / eps_ox * (tox - x) / tox
```

含义：

1. 每个缺陷族按一阶时间响应逐渐俘获电荷。
2. 每个缺陷族按自己的空间深度 `x_k` 屏蔽注入边界场强。
3. 注入边界场强 `E_inj(t)` 降低后，通过 FN 指数项使电流下降。

缺陷深度对电流的影响采用前面讨论的固定氧化层电压降近似：

```text
screen factor_k = (tox - x_k) / tox
E_inj(t) = E0 - sum_k[Q_def,k(t) / eps_ox * screen factor_k]
```

所以 `x` 越靠近注入界面，`screen factor` 越接近 1，电流衰减越强；`x` 越接近对侧界面，`screen factor` 越接近 0，对注入边界 FN 电流影响越弱。

当前 UI 中缺陷影响电流的链路为：

```text
x / E0 -> tau_eff -> Q_def(t) -> E_eff(t) -> I(t)
```

FN 系数和有效质量参数分工：

| 参数 | 当前作用 |
|---|---|
| `A_FN` | 常用 FN 公式中的前因子，`A_FN = q^3 m_SiC / (8π h m_ox Phi_B)` |
| `B_FN` | 指数衰减系数，`B_FN = 8π sqrt(2 m_ox Phi_B^3)/(3 q h)` |
| `Phi_B` | 固定为 SiC/SiO2 导带差 `2.7 eV` |
| `m_ox / m0` | 固定为 `0.29`，同时影响 `A_FN` 和 `B_FN`，其中对 `B_FN` 的指数影响最强 |
| `m_SiC / m0` | 固定为 `0.42`，进入 `A_FN` |
| `Area` | 将 `J_FN` 换算为总电流 `I = J_FN * Area` |

当前 UI 只使用物理 `A_FN` 公式计算电流密度，并通过 `Area` 换算为总电流；不再提供拟合 `K_FN` 前因子。

自动拟合实际识别的是等效缺陷族的电荷容量 `Q_k`、各缺陷族空间深度 `x_k` 和应力电压 `Vg`；`Phi_B` 固定为 `2.7 eV`，`mOx` 固定为 `0.29`，`mSiC` 固定为 `0.42`，`tau_eff` 由深度和局域场模型计算，不再把 `tau` 作为自由拟合参数。UI 中的 `Et` 固定为 `Ec(SiC)=0 eV`，不是从 Comphy 或 microscopic TAT 反演得到的唯一物理能级。

## 5. 时间常数模型

当前版本采用**空间深度与局域场耦合**的有效俘获时间常数模型。缺陷能级统一固定为 `Et=Ec(SiC)=0 eV`，不再进入 `tau_eff`，占据上限取 `f_inf=1`。因此三类缺陷的差异来自 `Nsheet`、输入 `tau` 和深度 `x`，不是来自不同能级。

基于表面电子浓度和隧穿深度的 reduced-order 计算流程为：

1. **场强与温度相关的动态表面势** <i>&psi;</i><sub>s</sub>(<i>E</i><sub>ox</sub>, <i>T</i>)：
   程序在正栅压积累区通过 Newton-Raphson 迭代求解泊松-玻尔兹曼方程，动态获取表面势 <i>&psi;</i><sub>s</sub>：
   ```text
   &psi;<sub>s</sub> = (k<sub>B</sub>T/q) * ln( (&epsilon;<sub>ox</sub><sup>2</sup> * E<sub>ox</sub><sup>2</sup>) / (2 * &epsilon;<sub>SiC</sub> * k<sub>B</sub>T * N<sub>D</sub>) )
   ```
   进而确定表面处费米能级位置与导带底的约化能量差 <i>&eta;</i><sub>s</sub> = (<i>E</i><sub>Fs</sub> - <i>E</i><sub>c</sub>(0)) / <i>k</i><sub>B</sub><i>T</i>。
2. **表面电子浓度** <i>n</i><sub>s</sub>(<i>E</i><sub>ox</sub>, <i>T</i>)：
   利用 1/2 阶费米-狄拉克积分计算简并与非简并条件下的表面电子浓度：
   ```text
   n<sub>s</sub> = N<sub>c</sub> * F<sub>1/2</sub>(&eta;<sub>s</sub>)
   ```
3. **有效时间常数** <i>&tau;</i><sub>eff</sub>：
   纯俘获时间常数随表面电子浓度 <i>n</i><sub>s</sub> 缩放，并用缺陷深度的指数项表示隧穿距离增加导致的俘获变慢：
   ```text
   &tau;<sub>eff</sub> = &tau;<sub>0</sub> * ( n<sub>s</sub>(Ref, 300 K) / n<sub>s</sub>(E<sub>ox</sub>, T) ) * exp(2 * gamma * x)
   ```
   其中 `x` 是缺陷到注入界面的深度，`gamma` 是左侧面板中的深度系数。当前实现不使用 `Et`、`Et_local` 或费米占据上限修正 `tau_eff`。
4. **缺陷电荷占据演变** <i>Q</i><sub>def</sub>(<i>t</i>)：
   当前占据上限固定为 1，每个缺陷族按一阶过程充电：
   ```text
   Q<sub>def</sub>(t) = Q<sub>cap</sub> * ( 1 - exp(-t / &tau;<sub>eff</sub>) )
   ```
这个处理保留了“缺陷越深，俘获越慢；场强和温度改变表面电子供给”的趋势，但没有引入能级相关的俘获/发射详细平衡。

为避免极端参数导致数值溢出，UI 对指数项和有效时间常数做了裁剪：

```text
-50 <= exponent <= 50
1e-12 s <= tau_eff <= 1e20 s
```

## 6. 局域场与势垒图

势垒图只画两条主曲线：

| 曲线 / 标记 | 含义 |
|---|---|
| 蓝色 `SiO2 导带 Uc(x)` | 氧化层导带 / 隧穿势垒高度，左轴单位 eV |
| 绿色 `氧化物场强 Eox(x)` | 局域氧化物场强，右轴单位 MV/cm |
| G1-G3 缺陷点 | 横坐标为缺陷深度 `x`，纵坐标固定为 `Et=Ec(SiC)=0 eV` |

当前版本的三类有效缺陷均固定在 `Et=Ec(SiC)=0 eV`。三类缺陷分别使用各自的 `Nsheet`、`tau` 和深度 `x`；`Et` 输入框为固定显示，不作为可调参数。

未俘获电荷时，可近似为均匀场三角势垒：

```text
U0(x) = Phi_B - integral_0^x E0 ds
```

存在多个深度处的俘获电荷后，UI 按 G1-G3 的位置 `x_k` 从左到右叠加局域电荷包，并在固定氧化层电压降条件下归一化局域场：

```text
integral_0^tox Eox(x) dx = E0 * tox
```

势垒由场强积分得到：

```text
U(x) = Phi_B - integral_0^x Eox(s) ds
```

这对应“俘获负电荷降低注入边界附近场强，使 FN 隧穿势垒变宽”的展示逻辑。

势垒图横轴固定为 `0-5 nm`，不是完整 `tox`，也不是按 `Uc=-3 eV` 自动截断。这样不同参数下可以在同一横轴范围内直观看到隧穿宽度 `xt` 的变化。

缺陷点的显示规则：

- 圆圈横坐标：各缺陷族自己的深度 `x`。
- 圆圈纵坐标：固定为 `Ec(SiC)=0 eV`。
- 圆圈大小：等效面密度 `Nsheet`，按 `sqrt(Nsheet / max Nsheet)` 缩放。
- 圆圈颜色深浅：当前占据 / 俘获程度 `fill(t)`。

势垒曲线仍会随局域场弯曲；缺陷能级标记不随势垒曲线移动。若只估算氧化层导带势能弯曲，在 `Eox ~= 9 MV/cm` 时，`1 nm` 对应的能量变化量级约为：

```text
Delta U ~= -Eox * 1 nm ~= -0.9 eV
```

所以界面附近 1 nm 内的导带弯曲并不小，但当前版本不把这部分弯曲转化为缺陷 `Et` 或 `tau_eff` 的变化。

注意：势垒图为了让局域弯曲更可见，对局域电荷显示做了视觉放大；这只影响图形展示，不改变 `I(t)`、`Q_def(t)` 或 KMC 的计算。

## 7. 缺陷浓度单位

当前 UI 中 G1-G3 的 `Nsheet` 单位是 `1e12 cm^-2`，表示从拟合电荷容量换算得到的等效面密度：

```text
Nsheet = Qcap / q
```

这是二维等效量，不是体陷阱浓度。若要换算为体浓度，需要给定对应缺陷族在氧化层中的有效厚度 `Delta x`：

```text
Nvol ~= Nsheet / Delta x
```

其中 `Nvol` 的单位为 `cm^-3`。当前 UI 不把 `Nsheet` 自动解释为体陷阱浓度，避免引入未校准的空间厚度假设。

## 8. 氧化层厚度起伏 Monte Carlo

该模块用于正向仿真不同器件间的氧化层微小起伏，而不是数据拟合。模型把每个器件拆成若干局域区域：

```text
tox_i = tox_mean + Delta tox_device + delta tox_local,i
E_i = Vg / tox_i
I_device(t) = sum_i J_FN(E_i, t) * Area_i
```

其中：

| 参数 | 含义 |
|---|---|
| `器件数` | Monte Carlo 中抽样的器件数量 |
| `每器件局域区域数` | 每个器件内部用于面积积分的局域 patch 数 |
| `σ_device` | 不同器件之间的整体厚度偏移标准差 |
| `σ_local` | 同一器件内部局域厚度起伏 RMS |
| `seed` | 固定随机序列，便于复现实验 |

主电流图中，紫色虚线为器件间中位数 P50，紫色阴影为 P5-P95。由于 FN 电流对局域场强呈指数敏感，较薄的局域 patch 会对总电流产生非线性放大，因此该模块直接对局域 `tox_i` 计算电流，而不是引入经验 `K_FN` 缩放因子。

同一批虚拟器件还会生成多条单器件 I-t 曲线，并根据器件平均厚度得到每个器件的 TDDB/KMC 击穿时间估计：

```text
E_device = Vg / mean(tox_i)
tBD_device = tBD_nominal * exp[-gamma_E * (E_device - E_nominal)]
```

其中 `tBD_nominal` 来自当前 KMC/TDDB E-model 的名义贯通时间。然后按排序失效概率：

```text
F_i = (i - 0.3) / (N + 0.4)
Y_i = ln[-ln(1 - F_i)]
X_i = ln(tBD_i)
```

绘制 Weibull 图，并给出线性拟合得到的 shape `beta`、characteristic life `eta` 和 `R^2`。这仍然是正向展示模型，不是对实验 Weibull 数据的拟合。

## 9. KMC 导电通道击穿演示

KMC 模块使用二维栅格表示氧化层截面，横向对应器件横向长度，纵向对应氧化层厚度：

```text
cell size = 1 nm x 1 nm
cols = round(device_length_nm / 1 nm)
rows = round(tox_nm / 1 nm)
```

因此，当 `tox = 50 nm` 时，KMC 深度方向使用约 50 行；横向行列数由左侧 `器件横向长度 (nm)` 控件给定。画布只是按比例缩放这个二维截面，不改变每个 KMC 单元代表的物理尺寸。

单元状态：

| 状态 | 图中颜色 | 含义 |
|---|---|---|
| 完整栅格 | 浅灰 | 尚未成为缺陷 |
| 初始缺陷 | 绿色 | 仿真开始前已有的随机缺陷 |
| KMC 生成缺陷 | 橙色 | 随 KMC 事件生成的新缺陷 |
| 贯通路径 | 红色 | gate 到 substrate 的连通通道 |

事件率采用 TDDB E-model：

```text
r_i(t) = k0 * exp[gamma_E * E_local,i(t)] * (1 + neighborBoost)^(n_i)
```

其中：

| 参数 | 含义 |
|---|---|
| `k0` | 基础缺陷生成率 |
| `E_local,i(t)` | 第 `i` 个 KMC 单元的局域氧化层场强，单位 MV/cm |
| `gamma_E` | E-model 场加速系数，单位 `(MV/cm)^-1` |
| `neighborBoost` | 相邻缺陷增强因子 |
| `n_i` | 当前单元相邻已缺陷单元数 |

`E_local,i(t)` 由当前 FN-defect 模型的注入边界有效场与各缺陷族的局域俘获电荷扰动近似得到。导入实验数据后，如果能识别到电流突跳击穿点，自动拟合会用实验 `tBD` 调整 `log10(k0)`，使 KMC 贯通时间接近实验击穿时间。`gamma_E` 和 `neighborBoost` 保持为用户可调参数，避免在 1 nm 物理网格上做过重的多参数 KMC 搜索。拟合完成后，页面会把当前时间同步到模型 `tBD`，因此主电流图中的击穿标记和 KMC 二维导电通道图展示的是同一个时间点。

KMC 时间推进：

```text
Delta t = -ln(u) / sum_i r_i
```

当缺陷从 gate 侧连通到 substrate 侧时，记录为：

```text
tBD = breakdown time
```

主电流图中的红色曲线是击穿后展示电流，用于直观标记 KMC 贯通后的跳变；当前版本校准的是 `tBD`，不是完整 post-breakdown 电导演化。

## 10. 自动保存与恢复参数

为了提高调试参数时的效率，页面引入了基于浏览器 `localStorage` 的参数自动暂存机制：
1. **自动保存**：每当您在左侧控制面板中输入新数值、勾选开关或拖动滑块时，修改后的最新参数集合都会实时序列化并以 `defect_tunneling_ui_params` 为键自动写入浏览器本地存储中。
2. **自动载入**：当您刷新页面、重新启动服务或更换端口再次进入时，UI 会在完成 DOM 初始化后自动读取该存储，并将修改后的参数自动还原填入控制面板，避免了每次打开都需要重新输入参数的繁琐步骤。

## 11. 已解决的核心物理模型更新 (能级相关俘获与高场能带弯曲)

为了使该可视化仿真的俘获时间常数模型与经典文献（如 Tewksbury & Lee 等）保持高度的物理一致性，当前版本已完成以下两点核心物理特性的升级：
1. **能级相关的多声子辅助晶格弛豫 (Lattice Relaxation) 修正**：
   * 当缺陷能级 <i>E</i><sub>t</sub> &lt; 0 eV 时，有效俘获截面引入晶格形变激活能 <i>E</i><sub>A</sub> 修正：<i>&sigma;</i><sub>eff</sub> = <i>&sigma;</i> &middot; exp(-<i>E</i><sub>A</sub> / <i>k</i>T)，其中 <i>E</i><sub>A</sub> = 0.15 &middot; |<i>E</i><sub>t</sub>| eV。这使得深能级缺陷的俘获速率相较浅能级指数级变慢，体现了非弹性隧穿的多声子阻碍。
2. **高场下能带弯曲导致的隧穿势垒降低 (Barrier Reduction)**：
   * 在 WKB 隧穿指数 &kappa; 计算中，考虑了电场对氧化层势垒的线性倾斜（梯形势垒平均高度降低项）。平均势垒高度修正为 <i>&Phi;</i><sub>avg</sub> = <i>&Phi;</i><sub>B</sub> - <i>E</i><sub>t</sub> - 0.05 &middot; <i>E</i><sub>eff</sub>(MV/cm) &middot; <i>x</i>(nm) eV。这真实地反映了高电场下隧穿概率随场强增大而指数级增加的物理效应。

## 12. 未完成的差异与后续扩展建议

若后续需要进一步提高拟合与物理仿真的精度，可着重考虑引入以下未完成的物理差异：

| 差异项 | 物理影响 | 扩展建议 |
|---|---|---|
| **空穴隧穿与双载流子参与** | 在极低电场或负 gate 偏置下，价带空穴隧穿机制主导，其势垒高度更大且有效质量 <i>m</i><sub>h</sub><sup>*</sup> 远大于电子。 | 当缺陷能级极深且场强反向时，切换为使用空穴参数（<i>m</i><sub>h</sub><sup>*</sup> &approx; 10 &middot; <i>m</i><sub>0</sub>）计算俘获时间常数。 |
| **界面态辅助两步过程 (&tau;<sub>it</sub>) 路径切换** | 电子在中性区或禁带中央的能级俘获，受限于热跃迁步骤（SRH 复合）。 | 引入能级相关的路径分段函数，在强积累区外自动将时间常数修正为 &tau;<sub>it</sub> = &tau;<sub>tun</sub> + &tau;<sub>SRH</sub>。 |
| **空间连续的缺陷能量与位置积分** | 当前采用 3-4 个离散“缺陷族”进行近似计算，无法完全匹配真实的体陷阱连续空间能量分布。 | 支持更细密的空间格点剖分（如 1D Poisson 求解格点），用沿 <i>x</i> 轴连续积分的电荷剖面代替离散电荷面。 |
| **自洽的能级空间移位 (Dynamic Band Bending)** | 局域电场 <i>E</i><sub>ox</sub>(<i>x</i>) 随时间被屏蔽，导致局域缺陷能级 <i>E</i><sub>t</sub> - <i>qE</i><sub>ox</sub><i>x</i> 随时间发生动态移动。 | 将时间常数 &tau;<sub>eff</sub> 放在时域递推循环中进行实时自洽更新（即使 &tau;<sub>eff</sub>(<i>t</i>) 成为时间的函数）。 |

## 13. 导出

页面支持：
- 导出当前电流图 PNG。
- 导出 CSV 曲线数据。

CSV 包含时间、电流、无缺陷基准、电流比值、`Einj`、`Qdef`、注入边界屏蔽因子、`Area`、`mOx`、`mSiC`、`A_FN`、`B_FN`、KMC 状态、KMC 击穿时间、KMC 二维截面尺寸、KMC 单元尺寸以及各缺陷族的 `x` 和 `tau_eff`。
