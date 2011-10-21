#!/usr/bin/env python
"""
Program to produce graphs from the output of the as_toroidal program
"""
from __future__ import division
from matplotlib.patches import Patch
from pylab import *
import sys
import math

# First argument to program is folder with results
folder=sys.argv[1]

# Open parameters file
params=open(folder+'params.dat','r')
data=params.readline().split(":");
rNum=int(data[1])
data=params.readline().split(":");
thNum=int(data[1])
data=params.readline().split(":");
dt=float(data[1])
data=params.readline().split(":");
tNum=int(data[1])
data=params.readline().split(":");
plotSteps=int(data[1])
data=params.readline().split(":");
rmin=float(data[1])
data=params.readline().split(":");
thtd=float(data[1])
params.close()

#solve values of dr and dth
dr=(1.0-rmin)/rNum
dth=math.pi/thNum

#Create array to store B at each step
B=zeros((rNum,thNum));

# make this smaller to increase the resolution
dx=0.005
# these variables are used by matplotlib
x = arange(0, 1, dx)
y = arange(-1, 1, dx)
X,Y = meshgrid(x, y)
Z=zeros((len(x),len(y)))

#get place in B grid that corresponds to each place in Z grid
grid=zeros((len(x),len(y),2))
i,j=0,0
for xvalue in x :
    for yvalue in y:
        r=math.sqrt(xvalue**2+yvalue**2)
        th=math.acos(yvalue/r)
        rStep=int((r-rmin)/dr)
        thStep=int(th/dth)
        if rStep<=0 or rStep>=rNum or thStep==0 or thStep>=thNum:
            grid[i][j][0]=0
            grid[i][j][1]=0
            j+=1
            continue
        grid[i][j][0]=rStep
        grid[i][j][1]=thStep
        j+=1
    j=0
    i+=1

k=0
while 1:
    #read first file
    data
    try:
        data=open(folder+"data_"+str(k),'r')
    except IOError as e:
        break
    i,j=0,0
    for line in data:
        values=line.split(" ")
        for value in values:
            if j==thNum:
                break
            B[i][j]=float(value)
            j+=1
        j=0
        i+=1
    data.close()

    i,j=0,0
    for xvalue in x :
        for yvalue in y:
            Z[i][j]=B[grid[i][j][0]][grid[i][j][1]]
            j+=1
        j=0
        i+=1
    #add data
    figtext(0.1, 0.8, "Radial steps: "+str(rNum))
    figtext(0.1, 0.7, "Angular steps: "+str(thNum))
    figtext(0.1, 0.6, "Time step: "+str(dt))
    figtext(0.1, 0.5, "t_h/t_d: "+str(thtd))
    figtext(0.1, 0.4, "t: "+str(k*dt))
    #create plot
    imshow(Z.T,extent=[0,1,-1,1],origin="lower")
    colorbar()
    #add zeros to the number of the plot, so they are ordered appropately
    num_file=str(k)
    diff_zeros=len(str(tNum))-len(str(k))
    while diff_zeros>0:
        num_file="0"+num_file
        diff_zeros-=1
    print num_file

    #save to file
    savefig(folder+"plot_"+num_file+".png", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1)
    close()
    k+=plotSteps

sys.exit()


