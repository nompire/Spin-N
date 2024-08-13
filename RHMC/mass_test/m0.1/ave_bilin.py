import sys
import glob
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

filename ="hmc_out"
plaq = []
for line in open(filename):
	if line.startswith('STOCH BILIN'):
		temp = line.split()
		plaq.append(float(temp[2]))
tot = np.sum(plaq)
N = len(plaq)
jk_plaq = np.empty(N,dtype=np.float64)
for i in range(N):
	temp = (tot- plaq[i])/(N-1.0)
	jk_plaq[i]=temp
ave_plaq = np.average(jk_plaq)
var_plaq = (N-1.0)*np.sum((jk_plaq-ave_plaq)**2) / float(N)
print(f'{ave_plaq:.4g} {var_plaq:.4g}')







