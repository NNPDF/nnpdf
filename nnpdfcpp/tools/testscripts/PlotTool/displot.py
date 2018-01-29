# nnpydf - nh 08/14

import matplotlib, sys, os, numpy
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from operator import add, sub
from math import log

# DIS style plot
class displot:
  def __init__(self,x,y, yerr, name, logx = False, logy = False, xlab = "x", leg = True):
    
    self.setname = name
    
    # Number of 'y' coordinates
    self.ny = len(x)
    self.xaxes = []
    self.yvals = []
    self.yerrs = []
    
    for ival in xrange(self.ny):
      self.xaxes.append(x[ival])
      self.yvals.append(y[ival])
      self.yerrs.append(yerr[ival])

    # Plot a legend?
    self.leg = leg
    
    self.fig = plt.figure()
    self.ax1 = self.fig.add_subplot('111')
    
    # Scales
    self.logx = logx
    if logx == True:
      self.ax1.set_xscale('log')
    self.logy = logy
    if logy == True:
      self.ax1.set_yscale('log')
    
    # Axis labels
    self.ax1.set_xlabel(xlab)
      
    # Gridlines
    self.ax1.xaxis.grid(True)
    self.ax1.yaxis.grid(True)

    # Plot title
    self.fig.text(0.13,0.92,name, fontsize=12)

  def addTheory(self, XV, CV, ER, name):
    
    inc = 0
    for i in xrange(self.ny):
      if len(CV[i]) != 0:
        yvalues = [x*(2**inc) for x in CV[i]]
        yerrors = [x*(2**inc) for x in ER[i]]
  
        # Higher and Lower boundaries
        CVup = map(add, yvalues, yerrors)
        CVdn = map(sub, yvalues, yerrors)
  
        if inc == 0:
          # plot CVs
          tmpplot1, = self.ax1.plot(XV[i], yvalues, label=name, linestyle='-')
          colour = tmpplot1.get_color()
        else:
          # plot CVs
          tmpplot1, = self.ax1.plot(XV[i], yvalues, label=None, linestyle='-', color = colour)
        
        # error bars
        self.ax1.fill_between(XV[i], CVdn, CVup, alpha=0.2,
                              facecolor=colour, linewidth=0)

        # increment positional counter 
        inc = inc + 1



  def savePlot(self,plotname):
    # Plot datapoints
    inc = 0
    for i in xrange(self.ny):
      if len(self.yvals[i]) != 0:
        yvalues = [x*(2**inc) for x in self.yvals[i]]
        yerrors = [x*(2**inc) for x in self.yerrs[i]]
        self.ax1.errorbar(self.xaxes[i], yvalues, yerr=yerrors,fmt='o', ecolor='black', linestyle='None', markersize=4, mfc='black')
        inc = inc + 1

    #Adjust axes
    xvals = sum(self.xaxes,[])
    if self.logx == False:
      xdiff = (5.0*(max(xvals) - min(xvals)))/100
    else:
      xdiff = 10.0*log(max(xvals) - min(xvals))/100
    
    axes=self.ax1.axis()
    
    if self.logx == False:
      self.ax1.axis([min(xvals) - xdiff, max(xvals) +xdiff, axes[2], axes[3]])
    #else:
    #  self.ax1.axis([min(xvals)*xdiff, max(xvals)/xdiff, axes[2], axes[3]])


    if self.leg == True:
      self.legend = self.ax1.legend(loc='best')
      self.legend.get_frame().set_alpha(0.8)
    
    self.fig.savefig(plotname+'.pdf')
    self.fig.savefig(plotname+'.png', dpi=80)
