# -*- coding: utf-8 -*-
import re

with open('defect_tunneling_ui.html', 'r', encoding='utf-8') as f:
    code = f.read()

code = code.replace('Qdef', 'sigmaTrap')
code = code.replace('Q_def', 'sigmaTrap')
code = code.replace('Q_k', 'sigmaTrap_k')
code = code.replace('Q (C)', 'sigmaTrap (C/cm2)')
code = code.replace('Q_total', 'sigmaTrap_total')
code = code.replace('Q_k *', 'sigmaTrap_k *')

code = code.replace('class="card" id="kmc_panel"', 'class="card" id="kmc_panel" style="display:none;"')
code = code.replace('class="card" id="weibull_panel"', 'class="card" id="weibull_panel" style="display:none;"')

with open('defect_tunneling_ui.html', 'w', encoding='utf-8') as f:
    f.write(code)
