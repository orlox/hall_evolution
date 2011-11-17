#!/usr/bin/env python
"""
Program to produce graphs from the output of the as_toroidal program
"""
import matplotlib.pyplot as plt
from pylab import *
import sys
import os
import math

#get folder name, if not specified as argument, use the latest one
folder=""
if len(sys.argv)<2:
    timestamp="0"
    for root, dirs, files in os.walk("."):
        for name in dirs:
            if name[:8]=="results_" and len(name)==18:
                newtimestamp=name[8:]
                if newtimestamp>timestamp:
                    folder=name
                    timestamp=newtimestamp
else:
    # First argument to program is folder with results
    folder=sys.argv[1]
#add "/" if it is not given in folder name
if folder[len(folder)-1]!="/":
    folder=folder+"/"

#Show what folder is being used
print("Plotting data in "+folder)

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
dr=(1.0-rmin)/(rNum-1)
dth=math.pi/(thNum-1)

#Create array to store A and B at each step
A=zeros((2*rNum-1,2*thNum-1));
B=zeros((rNum,thNum));

# make this smaller to increase the resolution
dx=0.005
# these variables are used by matplotlib
x = arange(0, 1, dx)
y = arange(-1, 1, dx)
X,Y = meshgrid(x, y)
ZA=zeros((len(x),len(y)))
ZB=zeros((len(x),len(y)))

#get place in A grid that corresponds to each place in Z grid
gridA=zeros((len(x),len(y),2))
i,j=0,0
for xvalue in x :
    for yvalue in y:
        r=math.sqrt(xvalue**2+yvalue**2)
        th=math.acos(yvalue/r)
        rStep=int((r-rmin+dr/4)/dr*2)
        thStep=int((th+dth/4)/dth*2)
        if rStep<=0 or rStep>=2*rNum-1 or thStep==0 or thStep>=2*thNum-1:
            gridA[i][j][0]=0
            gridA[i][j][1]=0
            j+=1
            continue
        gridA[i][j][0]=rStep
        gridA[i][j][1]=thStep
        j+=1
    j=0
    i+=1

#get place in B grid that corresponds to each place in Z grid
gridB=zeros((len(x),len(y),2))
i,j=0,0
for xvalue in x :
    for yvalue in y:
        r=math.sqrt(xvalue**2+yvalue**2)
        th=math.acos(yvalue/r)
        rStep=int((r-rmin+dr/2)/dr)
        thStep=int((th+dth/2)/dth)
        if rStep<=0 or rStep>=rNum or thStep==0 or thStep>=thNum:
            gridB[i][j][0]=0
            gridB[i][j][1]=0
            j+=1
            continue
        gridB[i][j][0]=rStep
        gridB[i][j][1]=thStep
        j+=1
    j=0
    i+=1

k=0
while 1:
    #do not plot this timestep if plot exists
    #add zeros to the number of the plot, so they are ordered appropately
    num_file=str(k)
    diff_zeros=len(str(tNum))-len(str(k))
    while diff_zeros>0:
        num_file="0"+num_file
        diff_zeros-=1
    try:
        data=open(folder+"plot_"+num_file+".png")
        data.close()
        k+=plotSteps
        continue
    except:
        pass

    #read A file
    data
    try:
        data=open(folder+"A_"+str(k),'r')
    except IOError as e:
        break
    i,j=0,0
    for line in data:
        values=line.split(" ")
        for value in values:
            if j==thNum:
                break
            A[i][j]=float(value)
            j+=1
        j=0
        i+=1
    data.close()
    #read B file
    try:
        data=open(folder+"B_"+str(k),'r')
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
            ZA[i][j]=A[gridA[i][j][0]][gridA[i][j][1]]
            ZB[i][j]=B[gridB[i][j][0]][gridB[i][j][1]]
            j+=1
        j=0
        i+=1
    fig=plt.figure()
    #add data
    figtext(0.05, 0.95, "Radial steps: "+str(rNum))
    figtext(0.3, 0.95, "Angular steps: "+str(thNum))
    figtext(0.55, 0.95, "Time step: "+str(dt))
    figtext(0.75, 0.95, "t_h/t_d: "+str(thtd))
    figtext(0.92, 0.95, "t: "+str(k*dt))
    #create plot
    a=fig.add_subplot(1,2,1)
    plt.imshow(ZA.T,extent=[0,1,-1,1],origin="lower")
    plt.colorbar()
    a.set_title("alpha")
    a=fig.add_subplot(1,2,2)
    plt.imshow(ZB.T,extent=[0,1,-1,1],origin="lower")
    plt.colorbar()
    a.set_title("beta")

    #save to file
    savefig(folder+"plot_"+num_file+".png", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1)
    close()
    print num_file

    k+=plotSteps

sys.exit()


