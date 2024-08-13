#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

L=4 
folders = 'm1*'
cfgs = []
cfgs2 = []
for folder_name in glob.glob(folders):
	cfg=str(folder_name)
	cfg=cfg.lstrip('m')
	if cfg not in cfgs:
		cfgs.append(float(cfg))
		cfgs2.append(cfg)
cfgs.sort()
cfgs2.sort()
print(cfgs)
print(cfgs2)
path = os.getcwd()

bin_size=20
a_plaq = []
plaq_err = []
outfile_name = 'ave_plaq_N'+str(L)+'.txt'
outfile = open(outfile_name,'w')


for Lambda in cfgs:	
	new_path = path+'/m'+str(Lambda)
	os.chdir(new_path)
	print(os.getcwd())
	filename ="hmc_out"
	plaq = []
	for line in open(filename):
		if line.startswith('STOCH BILIN'):
			temp = line.split()
			#plaq.append(np.sqrt(float(temp[6])*float(temp[6])+float(temp[7])*float(temp[7])))
			plaq.append(np.sqrt(float(temp[2])*float(temp[2])+float(temp[3])*float(temp[3])))
			#plaq.append(float(temp[2]))
	j=0
	bin_plaq = []
	tot_plaq=0
	for i in plaq:
		tot_plaq +=i
		if j==bin_size:
			bin_plaq.append(tot_plaq/bin_size)
			tot_plaq=0.0
			j=0
		j+=1
		
	tot = np.sum(bin_plaq)
	N = len(bin_plaq)
	jk_plaq = np.empty(N,dtype=np.float64)
	for i in range(N):
		temp = (tot- bin_plaq[i])/(N-1.0)
		jk_plaq[i]=temp
	ave_plaq = np.average(jk_plaq)
	var_plaq = (N-1.0)*np.sum((jk_plaq-ave_plaq)**2) / float(N)
	a_plaq.append(ave_plaq)
	plaq_err.append(var_plaq)
	os.chdir(path)
	print(f'{Lambda} {ave_plaq:.4g} {var_plaq:.4g}',file=outfile)
	print(f'{Lambda} {ave_plaq:.4g} {var_plaq:.4g}')

outfile.close()





#space =50 
fig,ax = plt.subplots()
ax.errorbar(cfgs,a_plaq, yerr=plaq_err,fmt='o',ms=5)
ax.set_ylabel(r'Plaquette')
ax.set_xlabel(r'$\lambda$')
#ax.xaxis.set_major_locator(ticker.MultipleLocator(space))

plt.show()  
sys.exit()
	
	
	
	




