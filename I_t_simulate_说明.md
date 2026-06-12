# SiC MOS 栅电流 Lucky-Defect 分析说明

## 1. 这份脚本现在做什么

`I_t_simulate.py` 现在不再使用原先那套“陷阱填充 + 有效场衰减 + 经验型 `K_FN / K_TAT`”拟合。

它已经改成和 `[lucky_defect_model.py](/E:/OneDrive/python_code/TDDB_Simulate/lucky_defect_model.py)` 一致的电流模型与参数：

- 相同的 SiC/SiO2 材料参数
- 相同的 FN 电流表达式
- 相同的两段 WKB lucky-defect `rho = J_TAT / J_FN`
- 相同的“一个器件由一个主导 lucky defect 控制局部增强电流”的思路

脚本有两种运行方式：

1. 传入 CSV，分析真实的预击穿电流轨迹
2. 不传参数，自动生成一条基于同一套 lucky-defect 公式的 demo 数据并直接跑通

---

## 2. 为什么要重写

旧版 `I_t_simulate.py` 的主要问题有三类：

1. 它不能直接运行，因为命令行强制要求 `csv_path`
2. 它的 FN/TAT 不是 `lucky_defect_model.py` 那套物理公式，而是另一套经验型时变模型
3. 说明文档里的公式虽然写了 LaTeX，但很多查看器对 `\[...\]` 渲染不稳定

所以现在做的改动是：

- CLI 允许无参数直接运行 demo
- 电流模型统一回 lucky-defect 主线
- 说明文档统一改成 `$$ ... $$` 数学块

---

## 3. 当前采用的物理模型

### 3.1 总电流分解

在恒定应力场下，当前脚本采用：

$$
J_{\mathrm{total}} = J_{\mathrm{FN}} + J_{\mathrm{TAT}}
$$

并且写成 lucky-defect 更常用的增强形式：

$$
J_{\mathrm{total}} = J_{\mathrm{FN}} \left( 1 + \rho \right)
$$

其中：

$$
\rho = \frac{J_{\mathrm{TAT}}}{J_{\mathrm{FN}}}
$$

这里的 `rho` 只由：

- 氧化层场 `E_ox`
- 势垒高度 `phi_b`
- 缺陷位置 `x`

决定。

当前版本不再引入：

- 陷阱占据率 `f(x,t)`
- `E_eff(t)` 场屏蔽
- `K_FN`
- `K_TAT`
- `eta * Qtrap / eps_ox`

---

## 4. 与 lucky_defect_model.py 保持一致的参数

脚本固定使用：

$$
\phi_b = 2.7\ \mathrm{eV}
$$

$$
m_{\mathrm{ox}} = 0.5\,m_0
$$

$$
m_{\mathrm{SiC}} = 0.39\,m_0
$$

这些数值和 `[lucky_defect_model.py](/E:/OneDrive/python_code/TDDB_Simulate/lucky_defect_model.py)` 一致。

---

## 5. FN 隧穿电流

### 5.1 FN 电流密度

FN 电流密度使用：

$$
J_{\mathrm{FN}} = A_{\mathrm{FN}} E_{\mathrm{ox}}^2
\exp\!\left(-\frac{B_{\mathrm{FN}}}{E_{\mathrm{ox}}}\right)
$$

其中

$$
A_{\mathrm{FN}}(\phi) =
\frac{q^3 m_{\mathrm{SiC}}}{16 \pi^2 \hbar m_{\mathrm{ox}} \phi}
$$

$$
B_{\mathrm{FN}}(\phi) =
\frac{4\sqrt{2m_{\mathrm{ox}}}}{3q\hbar}\,\phi^{3/2}
$$

脚本里真正用到的是 `phi = phi_b`，也就是：

$$
J_{\mathrm{FN}} =
A_{\mathrm{FN}}(\phi_b)\,E_{\mathrm{ox}}^2
\exp\!\left[
-\frac{1}{E_{\mathrm{ox}}}
\frac{4\sqrt{2m_{\mathrm{ox}}}}{3q\hbar}\phi_b^{3/2}
\right]
$$

代码函数：

- `fn_prefactor`
- `fn_action_coeff`
- `fn_current_density_a_per_cm2`

---

## 6. Lucky-defect 的两段 WKB

### 6.1 三角势垒宽度

总 FN 势垒宽度：

$$
W = \frac{\phi_b}{qE_{\mathrm{ox}}}
$$

代码函数：`barrier_width_nm`

---

### 6.2 sweet spot

2015 lucky-defect 模型中的 sweet spot 位置为：

$$
x_{\mathrm{sweet}} = W \left(1 - \frac{1}{\sqrt{2}}\right)
$$

代码函数：`sweet_spot_nm`

这和 `[lucky_defect_model.py](/E:/OneDrive/python_code/TDDB_Simulate/lucky_defect_model.py)` 里的写法一致。

---

### 6.3 局部势垒

对位于 `x` 处的缺陷，局部势垒高度为：

$$
\phi_x = \phi_b - qE_{\mathrm{ox}}x
$$

如果：

$$
x \le 0
$$

或

$$
x \ge W
$$

或

$$
\phi_x \le 0
$$

则当前缺陷对 lucky-defect TAT 不起作用，脚本直接返回：

$$
\rho = 0
$$

---

### 6.4 两段作用量

第一段从注入端到缺陷：

$$
S_1 = \frac{B}{E_{\mathrm{ox}}}
\left(\phi_b^{3/2} - \phi_x^{3/2}\right)
$$

第二段从缺陷到氧化层导带：

$$
S_2 = \frac{B}{E_{\mathrm{ox}}}\phi_x^{3/2}
$$

完整 defect-free FN 作用量：

$$
S_{\mathrm{FN}} = \frac{B}{E_{\mathrm{ox}}}\phi_b^{3/2}
$$

其中：

$$
B = \frac{4\sqrt{2m_{\mathrm{ox}}}}{3q\hbar}
$$

---

### 6.5 `rho = J_TAT / J_FN`

当前脚本按 `[lucky_defect_model.py](/E:/OneDrive/python_code/TDDB_Simulate/lucky_defect_model.py)` 的同一稳定写法计算：

$$
\log N = \log A_{\mathrm{FN}}(\phi_x) + (S_{\mathrm{FN}} - S_2)
$$

$$
\log D_1 = \log A_{\mathrm{FN}}(\phi_x) + (S_1 - S_2)
$$

$$
\log D_2 = \log A_{\mathrm{FN}}(\phi_b)
$$

$$
\log D = \log\!\left(e^{\log D_1} + e^{\log D_2}\right)
$$

$$
\rho = \exp(\log N - \log D)
$$

代码函数：`rho_tat_over_fn_nm`

这一步是本脚本和 `lucky_defect_model.py` 保持一致的核心。

---

## 7. 当前脚本怎样处理时间轴

这里必须明确：

- lucky-defect 原始模型不是“时变陷阱占据率模型”
- 它描述的是：在恒定应力下，一个主导 defect 会把 through-oxide current 提高

因此，当前脚本对真实 `I(t)` 数据采用的是“预击穿包络分解”方法，而不是完整的时变陷阱动力学拟合。

真实数据的处理顺序是：

1. 先确定击穿起点，只保留击穿前数据
2. 对击穿前电流做 rolling median 平滑，得到共同的预击穿包络
3. 用这个包络对应的代表性电流密度，反推出一个等效主导缺陷位置 `x`
4. 用 lucky-defect 模型得到固定的
   - `rho`
   - `alpha_FN = 1 / (1 + rho)`
   - `alpha_TAT = rho / (1 + rho)`
5. 再把同一条平滑包络按固定比例拆成
   - `I_FN(t) = alpha_FN \cdot I_env(t)`
   - `I_TAT(t) = alpha_TAT \cdot I_env(t)`

所以，这份脚本输出的总电流轨迹满足：

$$
I_{\mathrm{total}}(t) = I_{\mathrm{env}}(t)
$$

而分量分解满足：

$$
I_{\mathrm{FN}}(t) = \frac{1}{1+\rho} I_{\mathrm{env}}(t)
$$

$$
I_{\mathrm{TAT}}(t) = \frac{\rho}{1+\rho} I_{\mathrm{env}}(t)
$$

这一步保证了两件事同时成立：

- 电流分解所用的 `FN/TAT` 比例仍然来自 `lucky_defect_model.py`
- 模型总电流不会脱离真实预击穿包络

---

## 8. 击穿起点怎样确定

### 8.1 优先级

脚本按下面的优先级决定 `t_BD`：

1. 如果给了 `--breakdown-time-s`，直接使用这个时间
2. 如果输入是 Excel，且指定了 `--prefer-workbook-breakdown`，则优先读取 `Breakdown` sheet
3. 否则默认使用“第一电流急剧上升点”自动检测
4. 只有自动检测失败时，才回退到 workbook 里的 `BreakdownTime_s`

也就是说，当前默认行为已经是：

- **首个急剧上升点优先**

### 8.2 自动检测规则

自动检测函数是 `find_breakdown_index`。它扫描电流序列，找到第一个同时满足下面条件的点：

$$
I_i > I_{\mathrm{abs}}
$$

并且

$$
I_i > R_{\mathrm{BD}} \cdot \mathrm{median}\!\left(I_{i-window}, \dots, I_{i-1}\right)
$$

其中：

- `I_abs` 对应 `--breakdown-abs-a`
- `R_BD` 对应 `--breakdown-ratio`

默认值是：

- `breakdown_abs_a = 1e-5 A`
- `breakdown_ratio = 100`

对你这条 `CH04` 曲线，自动检测到的首个突升点正是：

$$
t_{BD} \approx 45397.220838\ \mathrm{s}
$$

对应电流约为：

$$
I_{BD} \approx 6.806643 \times 10^{-3}\ \mathrm{A}
$$

---

## 9. 等效缺陷位置怎么估计

### 9.1 代表性预击穿电流

脚本先对击穿前实测电流做 rolling median，得到：

$$
I_{\mathrm{env}}(t)
$$

再换算为电流密度：

$$
J_{\mathrm{env}}(t) = \frac{I_{\mathrm{env}}(t)}{A}
$$

然后取对数中位数作为代表值：

$$
J_{\mathrm{rep}} = 10^{\mathrm{median}\left(\log_{10} J_{\mathrm{env}}(t)\right)}
$$

### 9.2 扫描缺陷位置

接着在：

$$
0 < x < \min(W, t_{ox})
$$

范围内扫描位置，寻找使

$$
J_{\mathrm{ref}}(x) = J_{\mathrm{FN}} \left(1 + \rho(x)\right)
$$

最接近 `J_rep` 的 `x`。

如果：

$$
J_{\mathrm{rep}} \le J_{\mathrm{FN}}
$$

则脚本直接判为：

- `FN-only`
- 不引入 lucky defect
- `rho = 0`

这也是为什么当器件面积改为 `3 mm \times 3 mm` 后，单位面积电流密度明显下降，脚本会回到 `FN-only`。

---

## 10. 击穿前电荷积分

和旧版脚本一样，当前脚本仍然计算：

$$
Q(t_i) \approx \sum_{k=1}^{i}
\frac{I_{k-1}+I_k}{2}(t_k - t_{k-1})
$$

最终在最后一个预击穿点得到：

$$
Q_{BD} = \int_0^{t_{BD}^{-}} I(t)\,dt
$$

以及单位面积形式：

$$
\frac{Q_{BD}}{A} = \int_0^{t_{BD}^{-}} J(t)\,dt
$$

脚本同时输出四组：

1. 实测 `QBD`
2. 模型总电流 `QBD`
3. 纯 FN `QBD`
4. 纯 TAT `QBD`

---

## 11. 输入、输出与运行方式

### 11.1 输入

现在支持三类输入：

1. Keysight/B1500 风格 CSV

```text
DataValue, time, current, ...
```

2. 普通两列数值 CSV
3. Excel 工作簿 `xlsx/xls`

Excel 模式下可以用：

- `--sheet`
- `--time-column`
- `--current-column`

例如：

```bash
python I_t_simulate.py "PSU_D60_Run_20260401_164255.xlsx" ^
  --sheet I_t ^
  --time-column G ^
  --current-column H
```

### 11.2 输出

脚本会生成：

- `*_summary.txt`
- `*_fit_data.csv`
- `*_combined.png`
- `*_fn_only.png`
- `*_tat_only.png`
- `*_qbd.png`
- `*_position_scan.png`

其中 `summary.txt` 会明确写出：

- 击穿点来源是手工、自动首跳，还是 workbook
- 0 基索引和 1 基样本序号
- 原始 Excel/CSV 文件中的 1 基行号
- 预击穿包络的 rolling median 窗口
- `fit_status`

### 11.3 直接运行

无参数 demo：

```bash
python I_t_simulate.py
```

分析真实 Excel：

```bash
python I_t_simulate.py "PSU_D60_Run_20260401_164255.xlsx" ^
  --sheet I_t ^
  --time-column G ^
  --current-column H ^
  --vg 43 ^
  --tox-nm 45 ^
  --width-um 3000 ^
  --length-um 3000
```

如果已经知道击穿起点，直接固定：

```bash
python I_t_simulate.py "PSU_D60_Run_20260401_164255.xlsx" ^
  --sheet I_t ^
  --time-column G ^
  --current-column H ^
  --breakdown-time-s 45397.22084
```

如果已经知道要固定的缺陷位置，也可以直接指定：

```bash
python I_t_simulate.py "your_data.csv" --x-nm 0.35
```

---

## 12. 这份脚本的边界

它现在是：

- 与 lucky-defect-model 物理主线一致的电流分析脚本

但它不是：

- 完整的时变 TAT 俘获/释放仿真
- 完整的 `Et-x` 双变量能级分布仿真
- 2023 双缺陷带 Monte Carlo

更准确地说，这份脚本适合回答的问题是：

- 在给定场和氧化层厚度下，若预击穿电流主要由一个 lucky defect 主导，那么：
  - 纯 FN 电流有多大
  - TAT 增强电流有多大
  - 等效缺陷位置大概在哪
  - 到击穿前累计注入了多少电荷

---

## 13. 依赖

当前脚本只需要：

- `numpy`
- `matplotlib`

安装方式：

```bash
pip install numpy matplotlib
```
