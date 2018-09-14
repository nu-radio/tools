# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from pylab import setp

def beautify_limit(this_ax):

	"""
	beautify_limit

	Beautifies a limit plot


	Parameters
	----------
	this_ax (matplotlib.axes) : name of the axis you want beautified
		a matplotlib axis object

	Returns
	-------
	None:
		the function modifies the axes passed to it

	"""

	sizer=20
	xlow =  1.e15 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-20 #the lower y limit
	yup = 1e-10 #the upper y limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('E F(E) [$cm^{-2} s^{-1} sr^{-1}$]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_legend = this_ax.legend(loc='upper right',title='Single Event Sensitivity')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	this_ax.grid()

def beautify_aeff(this_ax):

	"""
	beautify_aeff

	Beautifies a effective area plot


	Parameters
	----------
	this_ax (matplotlib.axes) : name of the axis you want beautified
		a matplotlib axis object

	Returns
	-------
	None:
		the function modifies the axes passed to it

	"""

	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e2 #the lower x limit
	yup = 1.e10 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('[A$\Omega]_{eff}$  [cm$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='lower left',title='Effective Area')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)

def beautify_veff(this_ax):

	"""
	beautify_aeff

	Beautifies a effective volume plot


	Parameters
	----------
	this_ax (matplotlib.axes) : name of the axis you want beautified
		a matplotlib axis object

	Returns
	-------
	None:
		the function modifies the axes passed to it

	"""

	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e12 #the lower x limit
	yup = 1.e17 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('[V$\Omega]_{eff}$  [cm$^3$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='upper left',title='Effective Volume')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)

def beautify_counts(this_ax):

	"""
	beautify_counts

	Beautifies a histogram of the counts

	Parameters
	----------
	this_ax (matplotlib.axes) : name of the axis you want beautified
		a matplotlib axis object

	Returns
	-------
	None:
		the function modifies the axes passed to it

	"""

	sizer=20
	xlow =  1.e15 #the lower x limit
	xup = 1e21 #the uppper x limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('Events',size=sizer)
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='upper right',title='Event Counts')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)