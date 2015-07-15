from mpl_toolkits.mplot3d import axes3d
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import pdb,csv
from itertools import *
import brewer2mpl
import itertools
import numpy
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib import pyplot, mpl
from matplotlib.mlab import griddata
import math
from matplotlib import rcParams
from scipy.optimize import curve_fit
from scipy import polyval,linspace,polyfit

basePath = 'C:/Users/ritvik/Documents/PhD/Projects/CropIntensity/output/'
stateName = ['IA', 'ND', 'SD', 'MN', 'NE']
bmap = brewer2mpl.get_map('Set2', 'qualitative', 5, reverse=True)
colors = bmap.mpl_colors

index = 0
xMax = 0.0
yMax = 0.0
xMin = 0.0
yMin = 0.0

for st in stateName:
    
    fileName = basePath+'EXPR_FILE_'+st+'.csv'
    X = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(0))
    Y = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(1))
    Q = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(2))
    A1 = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(3))
    A2 = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(4))
    All = numpy.genfromtxt(fileName, delimiter=',',skip_header =1,usecols =(0,1,2,3,4))
    z=A1+A2
    
 
    # HEXAGON BINNING!!!
    plt.subplots_adjust(hspace=0.5)
    axes=plt.subplot(111)
    
    
    Q = Q+1.0
    Q = numpy.log10(Q)
    
    norm=mpl.colors.normalize(vmin=numpy.min(Q),vmax=numpy.max(Q))
    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'out'
    
                
    #plt.axis([z.min(), z.max()+1.0, 0.0, Q.max()])
    plt.title("Accuracy Vs Number of Rotations")
    #cb = plt.colorbar()
    #cb.set_label('counts')
    
    #rcParams['xtick.direction'] = 'out'
    #rcParams['ytick.direction'] = 'out'
    
    xlLim = z.min()
    xuLim = z.max()
    ylLim = Q.min()
    yuLim = Q.max()
    if(xlLim < xMin):
        xMin = xlLim
    if(xuLim > xMax):
        xMax = xuLim
    if(ylLim < yMin):
        yMin = ylLim
    if(yuLim > yMax):
        yMax = yuLim
            
    #axes.spines['right'].set_color('none')
    #axes.spines['top'].set_color('none')
    #axes.xaxis.set_ticks_position('bottom')
    

    
    # was: axes.spines['bottom'].set_position(('data',1.1*X.min()))
    #axes.spines['bottom'].set_position(('axes', -0.02))
    #axes.yaxis.set_ticks_position('left')
    #axes.spines['left'].set_position(('axes', -0.02))


    sc = plt.scatter(z, Q,marker='x',edgecolor='k',color=colors[index], s=60, alpha=0.65,lw=0.5)
    index+=1
    #etframes.add_dot_dash_plot(pylab.gca(), xs=z, ys=Q)
    #plt.colorbar(sc,orientation='vertical',shrink=1,extendfrac=0.1,aspect=50,fraction=0.15,spacing='uniform')

xrLim = (xMax-xMin)*0.15
yrLim = (yMax-yMin)*0.15
axes.set_xlim([xlLim-xrLim,xuLim+xrLim])
axes.set_ylim([ylLim-yrLim,yuLim+yrLim])
### GRID LINES!!!
axes.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75',alpha=0.5)
axes.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75',alpha=0.5)
axes.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75',alpha=0.5)
axes.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75',alpha=0.5)
            
plt.xlabel('Accuracy')
plt.ylabel('# of rotation (log)')
leg=plt.legend(('IA', 'ND', 'SD', 'MN', 'NE'),fancybox=False,loc='upper left')
leg.get_frame().set_alpha(0.5)
leg.draw_frame(False)
plt.savefig(basePath+'PlotALL'+'.png', bbox_inches=0)
