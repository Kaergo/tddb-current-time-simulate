# Current_simulate 说明

本目录用于 SiC MOS gate current 的 FN 隧穿、缺陷俘获、电流衰减、氧化层厚度起伏、KMC 导电通道击穿演示和预击穿 I-t 数据拟合。

当前可执行 UI：

```powershell
py -3.12 -S .\Current_simulate\run_defect_tunneling_ui.py
```

打开：

```text
http://127.0.0.1:8765/defect_tunneling_ui.html
```

## 1. 当前主模型

当前 UI 和脚本仍是 reduced-order 模型，不是 Comphy，也不是完整 microscopic TAT/NMP/Poisson 自洽求解器。

主链路为：

```text
defect charging -> Qdef(t) -> E_inj(t) -> FN current I(t)
```

电流模型：

```text
I(t) = Area * A_FN * E_inj(t)^2 * exp(-B_FN / E_inj(t))
```

缺陷俘获电荷通过屏蔽注入边界场强降低 FN 电流：

```text
\sigma_{trap}(t) = sum_k \sigma_{trap,k} * (1 - exp(-t / tau_eff,k))
E_inj(t) = E0 - sum_k[\sigma_{trap},k(t) / eps_ox * (tox - x_k) / tox]
```

其中 `x_k` 是缺陷到注入界面的深度。缺陷越靠近注入界面，`(tox - x_k) / tox` 越接近 1，对 FN 电流屏蔽越强。

固定物理参数建议：

```text
Phi_B = 2.7 eV
m_ox = 0.42 m0 或按当前实验设定固定
m_SiC = 0.29 m0 或按当前实验设定固定
```

## 2. 图表与 UI 功能

UI 主要包含：

```text
I(t) 电流曲线
tox Monte Carlo 多器件 I-t 与 Weibull 分布
KMC 2D 导电通道形成
SiO2 导带 Uc(x) 与氧化物局域场 Eox(x)
Qdef(t) 与 Eeff(t)
```

导入实验数据时，默认识别两列：

```text
CH11_Time_s
CH11_I_A
```

当前默认数据按 150 C 处理：

```text
T = 423.15 K
```

## 3. 预击穿拟合结果摘要

代表性输入：

```text
Current_simulate\1.xlsx
sheet = Sheet1
time column = CH11_Time_s
current column = CH11_I_A
tox = 42 nm
vg = 38 V
field ~= 9.0476 MV/cm
area = 4.0e-4 cm^2
```

自动检测到的击穿时间约为：

```text
tBD = 7.03053e5 s
```

旧版多时间常数等效拟合结果曾达到较高数值拟合精度：

```text
binned R2_log ~= 0.99994
raw R2_log ~= 0.99977
```

但这些 `\sigma_{trap,k}` 和 `tau_k` 只是 reduced-order 等效缺陷族参数，不能直接解释为唯一 microscopic trap species。

## 4. 截面与时间常数问题

早期 Scheme 3 曾尝试使用类似 NMP-inspired 的截面映射：

```text
sigma_eff(x) = sigma_N0 * exp(+/- 0.9 * E_depth(x) / kT)
```

该处理存在明显风险：

1. 若符号为正，截面可被指数放大到非物理量级，导致时间常数坍缩。
2. 若通过拟合直接压低截面，优化器可能得到极小的有效截面，例如 `~1e-33 cm^2`，这只是数学补偿，不是可信物理参数。
3. 在 `Eox ~= 9 MV/cm`、`sigma ~= 1e-15 cm^2`、`n ~= 1e19-1e20 cm^-3` 的量级下，近界面陷阱会非常快填充，难以解释长时间慢衰减。

因此后续不建议把 `sigma`、`Et`、`tau0` 和空间分布同时放开拟合。应固定或窄范围约束大部分物理量，只保留少数有明确物理意义的自由参数。

## 5. 能级与局域场的统一定义

后续能级模型应把缺陷本征能级定义为相对局部 SiO2 导带的深度：

```text
DeltaEt = Uc(x,t) - Etrap(x,t)
```

其中：

```text
Uc(x,t) = Phi_B - integral_0^x Eox(s,t) ds
Etrap_abs(x,t) = Uc(x,t) - DeltaEt
```

`DeltaEt` 对同一种氧化物缺陷固定不变；但 `Etrap_abs(x,t)` 会随氧化层局域场和能带弯曲改变。

正栅压下，电子势能随进入氧化层方向降低，因此远离界面时：

```text
Uc(x,t) 降低
Etrap_abs(x,t) 降低
```

不应把远离界面的缺陷能级画成相对 `Ec(SiC)` 反而升高。

## 6. 单能级拟合建议

当前若只拟合单温度、单一缺陷能级，建议采用：

```text
DeltaEt
Nt0
x_decay
```

作为主要拟合参数，并固定：

```text
Phi_B
m_ox
m_SiC
sigma 或 tau0
temperature
area
tox
```

空间分布：

```text
Nt(x) = Nt0 * exp(-x / x_decay)
```

时间推进中应每一步更新：

```text
Qtrap(x,t)
-> Eox(x,t)
-> Uc(x,t)
-> Etrap_abs(x,t)
-> tau_c(x,t)
-> I(t)
```

## 7. 双能级多温度拟合备选方案

后续如果有多温度 I-t 曲线，建议采用两个缺陷通道。

参考章节：

```text
Current_simulate\reference\1992\full1.md
  Chapter 3: Theory of Tunneling to and from Oxide Traps
  Section 3.4: Tunneling between Oxide Traps and Band States
  Section 3.4.1: Elastic Tunneling
  Section 3.4.2: Inelastic Tunneling
  Section 3.4.3: Lattice Relaxation Multiphonon Emission (LRME)
  Section 3.5: Tunneling between Oxide Traps and Interface States
  Section 4.3.1: The Effect of Band Bending in the Oxide
  Section 4.5: Temperature Dependence

Current_simulate\reference\1992\full2.md
  Appendix B: Tunneling Time Constants
  Appendix C: Hole Tunneling and Lattice Relaxation
  Appendix D: Solution of Coupled Rate Equations
  Appendix E: Band bending Calculations
```

导带对齐缺陷：

```text
tau_cb(x,t) = tau0_cb * exp[S_cb(x,t)]
```

禁带对齐缺陷：

```text
tau_inel(x,t,T) = tau0_inel(T) * exp[S_inel(x,t)]
```

导带对齐通道对应 SiC 导带电子弹性隧穿。禁带对齐通道不能简单等同为导带弹性隧穿，因为 SiC 禁带内没有同能量的连续电子态；参考 Tewksbury/Lee 模型，应理解为界面态辅助或晶格弛豫/多声子辅助的有效非弹性通道。

若界面态 SRH 交换很快：

```text
tau_it = tau_tun + tau_SRH ~= tau_tun
```

则该过程主要表现为不同的前因子，慢过程仍由隧穿指数控制。

与 1992 论文公式结构的对应关系：

```text
cb channel:
  Chapter 3.4.1 / Appendix B
  tau_cb = tau0_cb * exp(lambda_cb)

interface-trap assisted channel:
  Chapter 3.5 / Appendix D
  tau_it = tau_tun + tau_SRH
  当 tau_SRH 很短时，可近似为 tau_it ~= tau_tun

lattice-relaxation / multiphonon channel:
  Chapter 3.4.3 / Appendix C
  tau_lr = tau0_lr(T) * exp(lambda)
```

## 8. EA 的处理

`EA` 不是简单等于缺陷能级深度 `DeltaEt`。在晶格弛豫/NMP 图像中，`EA` 是构型坐标上的俘获激活能垒，受缺陷热力学能级、重组能、电子初态能量和局域电势共同影响。

参考章节：

```text
Current_simulate\reference\1992\full1.md
  Section 3.4.3: Lattice Relaxation Multiphonon Emission (LRME)
  Section 4.4.2: Lattice Relaxation
  Section 4.5.4: Lattice Relaxation

Current_simulate\reference\1992\full2.md
  Appendix C: Hole Tunneling and Lattice Relaxation
```

因此：

```text
EA != DeltaEt
EA = f(DeltaEt, relaxation energy, local electrostatic environment)
```

单温度下，`exp(EA/kT)` 与 `tau00_inel` 无法唯一分离，因此可以合并为常数前因子：

```text
tau0_inel_eff = tau00_inel * exp(EA/kT)
tau_inel = tau0_inel_eff * exp[S_inel]
```

多温度拟合时再显式引入：

```text
tau_inel(x,t,T) = tau00_inel * exp(EA/kT) * exp[S_inel(x,t)]
```

第一版建议把 `EA` 作为同一禁带对齐缺陷通道的常数参数，不随 `x` 和 `t` 动态变化。若多温度数据无法用常数 `EA` 解释，再考虑让 `EA` 依赖 `Etrap_abs(x,t)` 或局域场。

## 9. 晶格弛豫 (Lattice Relaxation) 公式与参数等效

根据以下核心文献，对于向 SiC 禁带深处（深陷阱）隧穿的过程，其物理图像和时间常数计算公式有明确的推导：
- [1] *Observation of Near-Interface Oxide Traps with the Charge-Pumping Technique* (R. E. Paulsen et al., IEEE TED 1992)
- [2] *Characterization, Modeling, and Minimization of Transient Threshold Voltage Shifts in MOSFET's* (T. L. Tewksbury et al.)

根据上述文献（包括 Transient Threshold Voltage Shifts 与 Charge Pumping 相关研究），向深能级隧穿的具体过程如下：

1. **两步走物理图像**：
   电子首先从 SiC 导带底（$E_c$）**水平地弹性隧穿**进入氧化层遇到缺陷的激发态；随后通过多声子发射（级联声子或晶格弛豫）**垂直向下**跌落至深基态。
   因此，WKB 隧穿积分 $\lambda$ 必须相对于导带底 $E_c$ 来计算（即代码中的 `Uc_local`），而不是相对于陷阱的绝对深度。

2. **多声子发射/晶格弛豫 (Lattice Relaxation, lr)**：
   $$ \tau_{lr} = \frac{1}{\sigma_0 \overline{v} [n_s + n_1 + p_s + p_1]} e^{E_A/kT} e^{\lambda(E_c, x)} $$
   在积累区 $n_s \gg n_1, p_s, p_1$，公式退化为：
   $$ \tau_{lr} = \left( \frac{e^{E_A/kT}}{\sigma_0 \overline{v} n_s} \right) \cdot e^{\lambda(E_c, x)} $$
   此过程不需要费米概率尾巴惩罚（$g_{access}$），电子直接从导带底越迁，但需付出热激活能垒 $E_A$ 的指数代价。

**当前的数学妥协（等效处理）**：
由于在单一应力与单一温度（如 150℃）下，$n_s$、$v_{th}$、未知的截面 $\sigma_0$ 以及热激活能 $E_A$ 均为全局常数。我们在模型中暂时将它们合并为一个单一的宏观独立前因子拟合参数 `tau0_inel`：
$$ \tau_{0,inel} \equiv \frac{e^{E_A/kT}}{\sigma_0 \overline{v} n_s} $$
代码实施为：`tau_c = tau0_inel / T_tunnel`。
此参数使得优化器可以自由调节前因子的量级，完美吸收深陷阱 $E_A$ 带来的指数级惩罚，从而解除了旧模型中导带电子概率对截面的刚性束缚。

> **【后续优化方向】**
> 目前将 $\sigma_0, n_s, E_A$ 全部打包进 `tau0_inel` 的做法物理上不够精细。后续在进行**多温度数据联合拟合**时，必须将该前因子解耦：
> 1. 显式引入实验温度 $T$。
> 2. 将 $E_A$ 作为独立的拟合参数，从而使得 $\exp(E_A/kT)$ 能随温度变化。
> 3. 显式代入 $n_s$ 和 $v_{th}$ 的具体数值，暴露出基础截面 $\sigma_0$ 进行拟合。

## 10. 多温度拟合参数建议

固定或窄范围约束：

```text
Phi_B = 2.7 eV
m_ox
m_SiC
sigma 或 tau0_cb
tox / area / stress voltage
```

优先拟合：

```text
DeltaEt_cb
DeltaEt_inel
Nt0_cb, x_decay_cb
Nt0_inel, x_decay_inel
tau0_cb 或窄范围 sigma_cb
tau00_inel
EA
```

如果参数仍然过多，可以先固定：

```text
x_decay_cb = x_decay_inel
```

或者固定 `tau0_cb`，优先识别温度依赖最强的 `EA`。

## 10. KMC 与厚度起伏模块

KMC 模块使用 TDDB E-model 演示缺陷生成与导电通道形成：

```text
r_i(t) = k0 * exp(gamma_E * E_local,i(t)) * (1 + neighborBoost)^n_i
```

二维网格按器件横向长度和氧化层厚度建立，默认单元为：

```text
1 nm x 1 nm
```

tox Monte Carlo 模块用于正向仿真不同器件间氧化层微小起伏导致的 FN 电流差异，并生成多条 I-t 曲线与 Weibull 分布。它不作为当前实验数据拟合的主要参数来源。

## 11. 相关文件

```text
defect_tunneling_ui.html                 可执行 UI
run_defect_tunneling_ui.py               本地 HTTP 启动脚本
I_t_simulate.py                          Python 拟合/仿真脚本
fit_prebreakdown.py                      预击穿拟合脚本
plot_tau_vs_x.py                         tau(x) 诊断图脚本
defect_tunneling_ui_说明.md              UI 详细说明
I_t_simulate_FN_defect_model说明.md      Python 模型历史说明
reference/                               参考文献与摘录
```
