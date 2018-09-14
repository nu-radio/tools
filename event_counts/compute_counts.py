# -*- coding: utf-8 -*-
import numpy as np
import sys
import matplotlib.pyplot as plt
from pylab import setp
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
from scipy.interpolate import splrep, splev

import constants as const
import tools as tool
import plots as plotter

'''
I want your energy in units of log10(eV)
And your veff in units of cm^3 * steradian
In a CSV file with one header line. See sample_veff.csv as an example.
'''

#stop them if they've forgotten a veff file
if(len(sys.argv)<2):
	print "Invalid usage! Use like: 'python compute_counts.py veff_file.csv'"
	print "Aborting!"
	exit()

filename = sys.argv[1]
data = np.genfromtxt(filename,delimiter=',',skip_header=1,names=['energy_logev','veff'])
logeV = data['energy_logev']
veff = data['veff']

#some more error handling for veff in the wrong units
if(np.max(veff)<1e10):
	print "Warning! Your veff*sr seem very small!"
	print "I suspect they are not in units of cm^3*sr"
	print "Please check this before trusting these results!"

#some more error handling
if(np.max(logeV)>25):
	print "Warning! Your energies seem very small!"
	print "I suspect they are not in units of log10(eV)"
	print "Please check this before trusting these results!"

aeff = veff/tool.get_Lint(np.power(10.,logeV))
interpolator = splrep(logeV, np.log10(aeff))

livetime = const.SecPerYear
counts_icecube=[]
counts_ahlers=[]
counts_koteramax=[]
energy_bins=[]
bins = np.arange(15.5,20.5,0.5)
for bin in bins:
	temp_logev = np.arange(bin,bin+0.5,0.1)
	temp_energy = np.power(10.,temp_logev)
	temp_aeff = np.power(10.,splev(temp_logev, interpolator))
	temp_icecube = tool.get_flux('icecube_thrumu',temp_logev)
	temp_ahlers = tool.get_flux('ahlers_2012',temp_logev)
	temp_koteramax = tool.get_flux('kotera_max',temp_logev)

	temp_counts_icecube = np.trapz(temp_icecube*temp_aeff*livetime,temp_energy)
	temp_counts_ahlers = np.trapz(temp_ahlers*temp_aeff*livetime,temp_energy)
	temp_counts_koteramax = np.trapz(temp_koteramax*temp_aeff*livetime,temp_energy)

	counts_icecube.append(temp_counts_icecube)
	counts_ahlers.append(temp_counts_ahlers)
	counts_koteramax.append(temp_counts_koteramax)
	energy_bins.append(np.power(10.,bin))

counts_icecube=np.array(counts_icecube)
counts_ahlers=np.array(counts_ahlers)
counts_koteramax=np.array(counts_koteramax)
energy_bins=np.array(energy_bins)

fig = plt.figure(figsize=(2.*11,8.5))
ax_veff = fig.add_subplot(1,2,1)
ax_counts = fig.add_subplot(1,2,2)
ax_veff.plot(np.power(10.,logeV),veff,'-o',linewidth=2.0,color='blue',label=r'ARIANNA SP Station')
n_icecube, bins, patches= ax_counts.hist(energy_bins,
									bins=np.power(10.,np.arange(15,22,0.5)),
									weights=counts_icecube,
									label=r'IceCube Thru-Mu E$^{-2.19}$: %.2f'%counts_icecube.sum(),
									fill=False, 
									stacked=True, 
									histtype='step', 
									edgecolor='blue',
									linewidth=4)
n_ahlers, bins, patches= ax_counts.hist(energy_bins,
									bins=np.power(10.,np.arange(15,22,0.5)),
									weights=counts_ahlers,
									label=r'Ahlers 2012: %.2f'%counts_ahlers.sum(),
									fill=False, 
									stacked=True, 
									histtype='step', 
									edgecolor='red',
									linewidth=4)
n_koteramax, bins, patches= ax_counts.hist(energy_bins,
									bins=np.power(10.,np.arange(15,22,0.5)),
									weights=counts_koteramax,
									label=r'Kotera Max: %.2f'%counts_koteramax.sum(),
									fill=False, 
									stacked=True, 
									histtype='step', 
									edgecolor='green',
									linewidth=4)

# print "Energy \t Counts IceCube \t Counts Ahlers "
# for this_bin, count_ic, count_ah in zip(bins,n_icecube,n_ahlers):
# 	print "%.2f \t %.4f \t %.4f"%(np.log10(this_bin), count_ic,count_ah)

plotter.beautify_veff(ax_veff)
plotter.beautify_counts(ax_counts)
ax_counts.set_ylabel('Events Per Year',size=17) #modify this axis title from plotter default
fig.savefig("example.png",edgecolor='none',bbox_inches="tight") #save the figure