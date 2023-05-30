#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:30:46 2023

@author: charmibhatt
"""

report = """
# fitting method   = leastsq
# function evals   = 567
# data points      = 300
# variables        = 39
chi-square         = 436.835953
reduced chi-square = 1.67370097
Akaike info crit   = 190.732577
Bayesian info crit = 335.180094
R-squared          = -9830.58357
[[Variables]]
B:         0.00248308 +/- 8.7988e-05 (3.54%) (init = 0.0023)
delta_B:  -0.06843322 +/- 0.00301885 (4.41%) (init = -0.0353)
zeta:     -0.31260631 +/- 0.00953055 (3.05%) (init = -0.4197)
T1:        84.4586063 +/- 6.82262452 (8.08%) (init = 67)
T2:        95.0389329 +/- 8.72110383 (9.18%) (init = 67)
T3:        98.2660464 +/- 9.76555277 (9.94%) (init = 67)
T4:        115.269670 +/- 22.1649526 (19.23%) (init = 67)
T5:        97.5633185 +/- 10.3067255 (10.56%) (init = 67)
T6:        84.7101410 +/- 9.68851068 (11.44%) (init = 67)
T7:        98.4130692 +/- 9.39173881 (9.54%) (init = 67)
T8:        88.6252830 +/- 6.69824023 (7.56%) (init = 67)
T9:        95.9090438 +/- 9.31651735 (9.71%) (init = 67)
T10:       87.5771631 +/- 7.14294904 (8.16%) (init = 67)
T11:       102.366810 +/- 10.5992689 (10.35%) (init = 67)
T12:       86.4965478 +/- 9.59656167 (11.09%) (init = 67)
sigma1:    0.18519313 +/- 0.00800038 (4.32%) (init = 0.2247)
sigma2:    0.20262523 +/- 0.00767339 (3.79%) (init = 0.2247)
sigma3:    0.19053018 +/- 0.00851895 (4.47%) (init = 0.2247)
sigma4:    0.19416707 +/- 0.01810849 (9.33%) (init = 0.2247)
sigma5:    0.21935156 +/- 0.00956987 (4.36%) (init = 0.2247)
sigma6:    0.16302652 +/- 0.01232314 (7.56%) (init = 0.2247)
sigma7:    0.17968789 +/- 0.00808328 (4.50%) (init = 0.2247)
sigma8:    0.19728615 +/- 0.00607026 (3.08%) (init = 0.2247)
sigma9:    0.20785736 +/- 0.00841500 (4.05%) (init = 0.2247)
sigma10:   0.22035833 +/- 0.00692019 (3.14%) (init = 0.2247)
sigma11:   0.24895832 +/- 0.00717047 (2.88%) (init = 0.2247)
sigma12:   0.19002270 +/- 0.01238269 (6.52%) (init = 0.2247)
origin1:   0.02905960 +/- 0.01003084 (34.52%) (init = 0.0041)
origin2:  -0.00949632 +/- 0.00755938 (79.60%) (init = 0.0041)
origin3:   0.00701112 +/- 0.00949471 (135.42%) (init = 0.0041)
origin4:  -0.00821518 +/- 0.02476964 (301.51%) (init = 0.0041)
origin5:   0.02816835 +/- 0.01115798 (39.61%) (init = 0.0041)
origin6:   1.8083e-04 +/- 0.01870953 (10346.55%) (init = 0.0041)
origin7:   0.06888491 +/- 0.00869012 (12.62%) (init = 0.0041)
origin8:   0.02449379 +/- 0.00447739 (18.28%) (init = 0.0041)
origin9:   0.12223918 +/- 0.00906592 (7.42%) (init = 0.0041)
origin10:  0.07789190 +/- 0.00600466 (7.71%) (init = 0.0041)
origin11:  0.04316670 +/- 0.00581988 (13.48%) (init = 0.0041)
origin12:  0.08014361 +/- 0.01787932 (22.31%) (init = 0.0041)
"""

# Save the report to a text file
# with open('fitting_12_sightlines.txt', 'w') as file:
#     file.write(report)
    



   # Define a class to store the variables
class Variables:
 def __init__(self):
     self.variables = {}

 def __getattr__(self, name):
     return self.variables.get(name, None)

# Create an instance of the Variables class
v = Variables()

# Open the lmfit report file and read its content
with open('fitting_12_sightlines.txt', 'r') as file:
    lines = file.readlines()
    parsing_variables = False
    for line in lines:
        if line.startswith("[[Variables]]"):
            parsing_variables = True
        elif parsing_variables and not line.startswith("#") and line.strip() != "":
            key_value = line.split(":")
            key = key_value[0].strip()
            if len(key_value) == 2:
                value = key_value[1].split("+/-")[0].strip()
                v.variables[key] = float(value)
            else:
                v.variables[key] = None

# Access the variables later
print(v.T1)  # Access T1
print(v.B)  # Access B
print(v.delta_B)  # Access delta_B
print(v.zeta)  # Access zeta
# ... access other variables as needed
