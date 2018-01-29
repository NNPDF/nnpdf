# nnpydf - nh 08/14

import matplotlib, sys, os, numpy, math
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from operator import add, sub
from math import log

# Kinematics ratio plot
class ratplot:
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
    
    # Scales
    self.logx = logx
    self.logy = logy

    # Setup gridspec
    gs = gridspec.GridSpec(self.ny, 1, height_ratios=[1 for x in self.xaxes])
    gs.update(wspace=0.00, hspace=0.00)
    
    self.fig = plt.figure()
    self.ax = []
    for g in gs:
      if len(self.ax) == 0:
        self.ax.append(self.fig.add_subplot(g))
      else:
        self.ax.append(self.fig.add_subplot(g, sharex=self.ax[0]))
      
      # Scales
      if logx == True:
        self.ax[-1].set_xscale('log')
      
      self.ax[-1].xaxis.grid(True)
      self.ax[-1].yaxis.grid(True)
        
      self.ax[-1].set_xlabel(xlab)
      
      plt.setp(self.ax[-1].get_xticklabels(), visible=False)


    # Plot title and ticks
    plt.setp(self.ax[-1].get_xticklabels(), visible=True)
    self.fig.text(0.13,0.92,name, fontsize=12)

  def addTheory(self, XV, CV, ER, name):
    
    inc = 0
    for i in xrange(self.ny):
      if len(CV[i]) != 0:
        yvalues = numpy.divide(CV[i],self.yvals[i])
        yerrors = numpy.divide(ER[i],self.yvals[i])

        # Higher and Lower boundaries
        CVup = map(add, yvalues, yerrors)
        CVdn = map(sub, yvalues, yerrors)
        
        if inc == 0:
          # plot CVs
          tmpplot1, = self.ax[i].plot(XV[i], yvalues, label=name, linestyle='-')
          colour = tmpplot1.get_color()
        else:
          # plot CVs
          tmpplot1, = self.ax[i].plot(XV[i], yvalues, label=None, linestyle='-', color = colour)
        
        # error bars
        self.ax[i].fill_between(XV[i], CVdn, CVup, alpha=0.2,
                              facecolor=colour, linewidth=0)
                              
        self.ax[i].yaxis.set_major_locator(MaxNLocator(4,prune='upper'))

        # increment positional counter 
        inc = inc + 1


  def savePlot(self,plotname):
    # Plot datapoints
    inc = 0
    for i in xrange(self.ny):
      if len(self.yvals[i]) != 0:
        yvalues = [1]*len(self.yvals[i])
        yerrors = numpy.divide(self.yerrs[i],self.yvals[i])
        
        self.ax[i].errorbar(self.xaxes[i], yvalues, yerr=yerrors,fmt='o', ecolor='black', linestyle='None', markersize=4, mfc='black')
        inc = inc + 1

    #if self.leg == True:
      #self.legend = self.ax[0].legend(loc='best')
      #self.legend.get_frame().set_alpha(0.8)
    
    #Adjust axes
    xvals = sum(self.xaxes,[])
    if self.logx == False:
      xdiff = (5.0*(max(xvals) - min(xvals)))/100
    else:
      xdiff = 10.0*log(max(xvals) - min(xvals))/100
    
    axes=self.ax[0].axis()
    
    if self.logx == False:
      self.ax[0].axis([min(xvals) - xdiff, max(xvals) +xdiff, axes[2], axes[3]])
    #else:
    #  self.ax[0].axis([min(xvals)*xdiff, max(xvals)/xdiff, axes[2], axes[3]])
    
    self.fig.savefig(plotname+'.pdf')
    self.fig.savefig(plotname+'.png', dpi=80)
