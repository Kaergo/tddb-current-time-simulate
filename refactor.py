import re
import os

filepath = 'I_t_simulate.py'
with open(filepath, 'r', encoding='utf-8') as f:
    code = f.read()

replacements = [
    (r'\bQ_def\b', 'sigma_trap'),
    (r'\bqdef\b', 'sigma_trap'),
    (r'\bQ_k\b', 'sigma_trap_k'),
    (r'\bqcaps\b', 'sigma_sat'),
    (r'\bqcaps_list\b', 'sigma_sat_list'),
    (r'\bQ_cap_max_x\b', 'sigma_sat_max_x'),
    (r'\bqeq\b', 'sigma_eq'),
    (r'\bQdef\b', 'sigmaTrap'),
]

for old, new in replacements:
    code = re.sub(old, new, code)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(code)

print('Basic renames done in I_t_simulate.py')
