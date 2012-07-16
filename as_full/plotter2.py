#!/usr/bin/env python
"""
Program to produce graphs from the output of the as_toroidal program
"""
import matplotlib.pyplot as plt
from pylab import *
import sys
import os
import math
from scipy import special
from matplotlib.colors import LinearSegmentedColormap

#Use a monospace font to avoid change in the size of plots
#rc('font',**{'family':'monospace','monospace':['Computer Modern Typewriter']})
rc('font',**{'family':'monospace'})

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
factor=float(data[1])
data=params.readline().split(":");
tNum=int(data[1])
data=params.readline().split(":");
plotSteps=int(data[1])
data=params.readline().split(":");
rmin=float(data[1])
data=params.readline().split(":");
thtd=float(data[1])
data=params.readline().split(":");
lNum=int(data[1])
params.close()

#solve values of dr and dth
dr=(1.0-rmin)/(rNum-1)
dth=math.pi/(thNum-1)

#Create array to store A and B at each step, additional steps on A are to plot field lines outside the star
A=zeros((2*rNum,thNum));
B=zeros((rNum,thNum));

#Create arrays that store cartesian coordinates of each point of the grid
XBvalues=zeros((rNum,thNum));
YBvalues=zeros((rNum,thNum));
for i in range(rNum):
    #use twice the value of the radial size of the crust to help visualization
    r=i*2*dr+(2*rmin-1)
    for j in range(thNum):
        th=j*dth
        XBvalues[i][j]=r*math.sin(th)
        YBvalues[i][j]=r*math.cos(th)
XAvalues=zeros((2*rNum,thNum));
YAvalues=zeros((2*rNum,thNum));
for i in range(2*rNum):
    #use twice the value of the radial size of the crust to help visualization
    r=0
    if i<=rNum-1:
        r=i*2*dr+(2*rmin-1)
    else:
        r=1+(i-rNum+1)*2.0/rNum
    for j in range(thNum):
        th=j*dth
        XAvalues[i][j]=r*math.sin(th)
        YAvalues[i][j]=r*math.cos(th)

#Create array for multipole coefficients
multipoles=zeros((lNum));
#store P_l^1(cos(theta)) at the surface, used to solve multipoles. I have an array evaluated in midpoints to perform
#the integrations, and another one solve on mesh points to perform alpha evaluation
plone=zeros((lNum,thNum));
plone2=zeros((lNum,thNum));
for j in range(thNum):
    alp=special.lpmn(1, lNum+1, math.cos((j+0.5)*dth))
    alp2=special.lpmn(1, lNum+1, math.cos(j*dth))
    for l in range(lNum):
        plone[l][j]=alp[0][1][l+1]
        plone2[l][j]=alp2[0][1][l+1]

#Create gray-scale printable colormap, cdict taken from user Damon McDougall-2 who posted it in the matplotlib mailing list
cdict = {'red':((0.000, 0.00, 0.00), 
                (0.125, 0.15, 0.15), 
                (0.250, 0.30, 0.30), 
                (0.375, 0.60, 0.60), 
                (0.500, 1.00, 1.00), 
                (0.625, 0.90, 0.90), 
                (0.750, 0.90, 0.90), 
                (1.000, 0.90, 0.90)), 
       'green':((0.000, 0.00, 0.00), 
                (0.125, 0.15, 0.15), 
                (0.250, 0.15, 0.15), 
                (0.375, 0.20, 0.20), 
                (0.500, 0.25, 0.25), 
                (0.625, 0.50, 0.50), 
                (0.750, 0.75, 0.75), 
                (1.000, 0.75, 0.75)), 
        'blue':((0.000, 0.00, 0.00), 
                (0.125, 0.50, 0.50), 
                (0.250, 0.75, 0.75), 
                (0.375, 0.50, 0.50), 
                (0.500, 0.15, 0.15), 
                (0.625, 0.00, 0.00), 
                (0.750, 0.10, 0.10), 
                (1.000, 0.10, 0.10))} 
mycmap=LinearSegmentedColormap('MyCpam', cdict)

k=0
maxA=0.0
maxB=0.0
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

    #read B file
    try:
        data=open(folder+"B_"+str(k),'r')
    except IOError as e:
        break
    #first line has simulation time
    t=data.readline()
    i,j=0,0
    for line in data:
        values=line.split(" ")
        for value in values:
            if j==thNum:
                break
            B[i][j]=float(value)
            if(math.fabs(B[i][j])>maxB):
                maxB=math.fabs(B[i][j])
            j+=1
        j=0
        i+=1
    data.close()

    #read A file
    data
    try:
        data=open(folder+"A_"+str(k),'r')
    except IOError as e:
        break
    #first line has simulation time
    t=data.readline()
    i,j=0,0
    for line in data:
        values=line.split(" ")
        for value in values:
            if j==thNum:
                break
            A[i][j]=float(value)
            if(math.fabs(A[i][j])>maxA):
                maxA=math.fabs(A[i][j])
            j+=1
        j=0
        i+=1
    data.close()
    #solve multipoles and evaluate alpha outside the star
    for l in range(lNum):
        multipoles[l]=0
        for j in range(0,thNum-1):
            multipoles[l]+=(A[rNum-1][j]+A[rNum-1][j+1])/2*plone[l][j]*dth
        multipoles[l]=multipoles[l]*(2*l+3)/2/(l+2)/(l+1)
    for i in range(rNum-1,2*rNum):
        r=1+(i-rNum+1)*1.0/rNum
        for j in range(0,thNum):
            th=j*dth
            A[i][j]=0
            for l in range(lNum):
                A[i][j]+=multipoles[l]/pow(r,l+1)*math.sin(th)*plone2[l][j]
    

    #set contour levels with the zero value always as the mid-contour
    numlevels=21
    levelsB=zeros([numlevels],float)
    levelsB[(numlevels-1)/2]=0
    levelsB[numlevels-1]=maxB
    levelsB[0]=-maxB
    for i in range(1,10):
        levelsB[(numlevels-1)/2+i]=maxB*i/(numlevels-1)*2
        levelsB[(numlevels-1)/2-i]=-maxB*i/(numlevels-1)*2
    levelsA=zeros([numlevels],float)
    levelsA[(numlevels-1)/2]=0
    levelsA[numlevels-1]=maxA
    levelsA[0]=-maxA
    for i in range(1,10):
        levelsA[(numlevels-1)/2+i]=maxA*i/(numlevels-1)*2
        levelsA[(numlevels-1)/2-i]=-maxA*i/(numlevels-1)*2

    fig=plt.figure()
    plt.ylim((-1,1))
    plt.xlim((0,2))
    #add data
    #figtext(0.02, 0.95, "Radial steps: "+str(rNum)+"     Angular steps: "+str(thNum)+"     Factor: "+str(factor)+ "     t_h/t_d: "+str(thtd)+"     t: "+t)
    figtext(0.35, 0.9, "t/t_h: "+t)
    #create plot
    plt.contourf(XBvalues,YBvalues,B,levels=levelsB, cmap=mycmap)
    plt.colorbar(format='%.2e')
    #Hide x and y axes
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)

    #plt.contour(XAvalues,YAvalues,A,levels=levelsA, colors="#33aa33")
    #Replacing the line just before one with the next two includes labels
    Aplot=plt.contour(XAvalues,YAvalues,A,levels=levelsA, colors="#33aa33")
    plt.clabel(Aplot, inline=1, fontsize=10, fmt="%.2e")

    #save to file
    savefig(folder+"plot_"+num_file+".png", dpi=None, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches='tight')
    close()
    print num_file

    #reset maximun values to zero so they are recalculated next step
    maxA=0
    maxB=0

    k+=plotSteps

sys.exit()


