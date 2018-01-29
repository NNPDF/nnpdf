# nnpydf - nh 08/14

import matplotlib, sys, os, numpy
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from operator import add, sub

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

# Standard Data/Theory comparison
class stdplot:
  def __init__(self,x,y, yerr, name, logx = False, logy = False, xlab = "x", leg = True):
    
    self.xaxis = list(x)
    self.data = list(y)
    self.error = list(yerr)

    # Plot a legend?
    self.leg = leg
    
    # Setup gridspec
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    gs.update(wspace=0.00, hspace=0.00)
    
    self.fig = Figure()
    self.canvas = FigureCanvas(self.fig)

    self.ax1 = self.fig.add_subplot(gs[0])
    self.ax2 = self.fig.add_subplot(gs[1])
    
    # Logarithmic x axis
    self.logx = logx
    if self.logx == True:
      self.ax1.set_xscale('log')
      self.ax2.set_xscale('log')

    # Logarithmic y axis
    self.logy = logy
    if self.logy == True:
      self.ax1.set_yscale('log')

    # Disable plot x ticks
    plt.setp(self.ax1.get_xticklabels(), visible=False)
    
    # Axis formatting
    self.ax2.yaxis.tick_right()
    self.ax2.yaxis.set_major_locator(MaxNLocator(5,prune='upper'))

    # Axis labels
    self.ax2.set_xlabel(xlab)

    # gridlines
    self.ax1.xaxis.grid(True)
    self.ax1.yaxis.grid(True)
    self.ax2.xaxis.grid(True)
    self.ax2.yaxis.grid(True)
      
    # Plot title
    self.fig.text(0.13,0.94,name, fontsize=12)

  def addTheory(self,CV, ER, name):
    # Higher and Lower boundaries
    CVup = map(add, CV, ER)
    CVdn = map(sub, CV, ER)
    
    
    if len(CV) > 1:
      # plot CVs
      tmpplot1, = self.ax1.plot(self.xaxis, CV, label=name)
      tmpplot2, = self.ax2.plot(self.xaxis, numpy.divide(CV,self.data), label=name)
      
      # error bars
      self.ax1.fill_between(self.xaxis, CVdn, CVup, alpha=0.2,
                            facecolor=tmpplot1.get_color(), linewidth=0)
                            
      # error bars
      self.ax2.fill_between(self.xaxis, numpy.divide(CVdn,self.data), numpy.divide(CVup,self.data),
                            alpha=0.2, facecolor=tmpplot2.get_color(), linewidth = 0)
    else:
      bar1 = self.ax1.errorbar(self.xaxis, CV, ER, fmt = 'x', markersize=4, label=name)
      self.ax2.errorbar(self.xaxis, numpy.divide(CV,self.data), numpy.divide(ER,self.data), fmt = 'x', markersize=4, label=name)


  def savePlot(self,plotname):
    # Plot datapoints
    self.ax1.errorbar(self.xaxis, self.data, yerr=self.error, fmt='o', ecolor='black', linestyle='None', markersize=4, mfc='black')
    self.ax2.errorbar(self.xaxis, [1]*len(self.data), yerr=numpy.divide(self.error,self.data), fmt='o', ecolor='black', linestyle='None', markersize=4,  mfc='black')
    
    # Adjust axes
    if len(self.data) > 2:
      xdiff = 5.0*((self.xaxis[-1] - self.xaxis[0])/100.0)
      axes=self.ax1.axis()
      self.ax1.axis([self.xaxis[0] - xdiff, self.xaxis[-1] +xdiff, axes[2], axes[3]])
      axes2=self.ax2.axis()
      self.ax2.axis([self.xaxis[0] - xdiff, self.xaxis[-1] +xdiff, axes2[2], axes2[3]])

    #if self.leg == True:
    #  self.legend = self.ax1.legend(loc='best')
    #  self.legend.get_frame().set_alpha(0.8)

    self.fig.savefig(plotname+'.pdf')
    self.fig.savefig(plotname+'.png', dpi=80)


  def close(self):
    self.fig.close()