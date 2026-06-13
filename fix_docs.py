# -*- coding: utf-8 -*-
import re
import os

files = ['README.md', 'I_t_simulate_FN_defect_model说明.md', 'I_t_simulate_说明.md']

for filepath in files:
    if not os.path.exists(filepath): continue
    with open(filepath, 'r', encoding='utf-8') as f:
        code = f.read()

    # Title replacements
    code = code.replace('SiC MOS 缺陷充电与 TDDB 电流模拟程序', 'FN Tunneling Current Relaxation Induced by Oxide-Trap Charging')
    code = code.replace('SiC MOS Defect Charging & TDDB Current Simulator', 'FN Tunneling Current Relaxation Induced by Oxide-Trap Charging')
    code = code.replace('SiC MOS 栅电流 I-t 数据模拟器', 'FN Tunneling Current Relaxation Simulator')

    # Variable replacements
    code = code.replace('Q_{def}', '\sigma_{trap}')
    code = code.replace('Q_def', '\sigma_{trap}')
    code = code.replace('Q_{def}(t)', '\sigma_{trap}(t)')
    code = code.replace('Q_k', '\sigma_{trap,k}')
    code = code.replace('Q_{total}(t)', '\sigma_{trap,total}(t)')
    
    # Equation replacements
    code = code.replace('\frac{Q_{def,k}(t)}{\varepsilon_{ox}}', '\frac{\sigma_{trap,k}(t)}{\varepsilon_{ox}}')
    code = code.replace('Q_def(t)', 'sigma_trap(t)')

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(code)

print('Documentation updated.')
