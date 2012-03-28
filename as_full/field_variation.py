#!/usr/bin/env python
"""
Program to compare how much the structure of a field changes with respect to its initial form
"""
import matplotlib.pyplot as plt
from pylab import *
import sys
import os
import math

#Fundamental Ohm mode for rmin=0.75
A1=-0.52004
A2=-0.55882
kk=7.03266
def A_Ohm_function(r,th):
    return r*pow(sin(th),2)*(((sin(kk*r)/(kk*r)-cos(kk*r))*A2)/(kk*r)+((-sin(kk*r)-cos(kk*r)/(kk*r))*A1)/(kk*r))

#Get name of folder with results
folder=sys.argv[1]
#add "/" if it is not given in folder name
if folder[len(folder)-1]!="/":
    folder=folder+"/"

#open output file and print header
f = open(folder+"compare.dat", 'w')
f.write("#num   t/t_h   int(dB)\n")

#Show what folder is being used
print("Analyzing data in "+folder)

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
A=zeros((rNum+1,thNum+1));
A_Ohm=zeros((rNum+1,thNum+1));
B=zeros((rNum+1,thNum+1));

#Fill array with ohm eigenmode for A_Ohm
for i in range(rNum):
    r=i*dr+rmin
    for j in range(thNum):
        th=j*dth
        A_Ohm[i][j]=A_Ohm_function(r,th)

#Create array to store the vector values of the field at the initial and later times
#field is solved at midpoints in the grid.
B_field_i=zeros((rNum,thNum,3));
B_field_k=zeros((rNum,thNum,3));
B_field_Ohm=zeros((rNum,thNum,3));
initial_energy=0;

#solve vector magnetic field of the fundamental Ohm mode
for i in range(rNum-1):
    r=i*dr+dr/2+rmin
    for j in range(thNum-1):
        th=j*dth+dth/2
        #Solve each component
        #r component
        B_field_Ohm[i][j][0]=1/r/r/sin(th)*(A_Ohm[i][j+1]-A_Ohm[i][j]+A_Ohm[i+1][j+1]-A_Ohm[i+1][j])/2/dth
        #th component
        B_field_Ohm[i][j][1]=-1/r/sin(th)*(A_Ohm[i+1][j]-A_Ohm[i][j]+A_Ohm[i+1][j+1]-A_Ohm[i][j+1])/2/dr
        #phi component
        B_field_Ohm[i][j][2]=0

#solve energy of Ohm eigenmode
energy_Ohm=0
for i in range(rNum-1):
    r=i*dr+dr/2+rmin
    for j in range(thNum-1):
        th=j*dth+dth/2
        energy_Ohm+=(pow(B_field_Ohm[i][j][0],2)+pow(B_field_Ohm[i][j][1],2)+pow(B_field_Ohm[i][j][2],2))*pow(r,2)*sin(th)
energy_Ohm=energy_Ohm*dr*dth/4

#analyze all data
k=0
while 1:
    #add zeros to the number of the plot, so they are ordered appropately
    num_file=str(k)
    diff_zeros=len(str(tNum))-len(str(k))
    while diff_zeros>0:
        num_file="0"+num_file
        diff_zeros-=1

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

    Bmaxth=0
    #solve vector magnetic field
    for i in range(rNum-1):
        r=i*dr+dr/2+rmin
        for j in range(thNum-1):
            th=j*dth+dth/2
            #Solve each component
            #r component
            B_field_k[i][j][0]=1/r/r/sin(th)*(A[i][j+1]-A[i][j]+A[i+1][j+1]-A[i+1][j])/2/dth
            #th component
            B_field_k[i][j][1]=-1/r/sin(th)*(A[i+1][j]-A[i][j]+A[i+1][j+1]-A[i][j+1])/2/dr
            #phi component
            B_field_k[i][j][2]=1/r/sin(th)*(B[i][j]+B[i+1][j]+B[i][j+1]+B[i+1][j+1])/4


    #solve total energy
    energy=0
    for i in range(rNum-1):
        r=i*dr+dr/2+rmin
        for j in range(thNum-1):
            th=j*dth+dth/2
            energy+=(pow(B_field_k[i][j][0],2)+pow(B_field_k[i][j][1],2)+pow(B_field_k[i][j][2],2))*pow(r,2)*sin(th)
    energy=energy*dr*dth/4

    #if this is the first timestep, store initial_energy and field.
    if k==0:
        initial_energy=energy
        for i in range(rNum-1):
            for j in range(thNum-1):
                B_field_i[i][j][0]=B_field_k[i][j][0]
                B_field_i[i][j][1]=B_field_k[i][j][1]
                B_field_i[i][j][2]=B_field_k[i][j][2]
    
    #solve integral of deltaB, at each time, field must be corrected to have equal energy than the initial one
    dB_energy=0
    for i in range(rNum-1):
        r=i*dr+dr/2+rmin
        for j in range(thNum-1):
            th=j*dth+dth/2
            dB_energy+=(pow(B_field_i[i][j][0]*sqrt(energy_Ohm/initial_energy)-B_field_k[i][j][0]*sqrt(energy_Ohm/energy),2)
                            +pow(B_field_i[i][j][1]*sqrt(energy_Ohm/initial_energy)-B_field_k[i][j][1]*sqrt(energy_Ohm/energy),2)
                            +pow(B_field_i[i][j][2]*sqrt(energy_Ohm/initial_energy)-B_field_k[i][j][2]*sqrt(energy_Ohm/energy),2))*pow(r,2)*sin(th)
    dB_energy=dB_energy*dth*dr/4

    #solve integral of deltaB with respect to Ohm fundamental mode, at each time, field must be corrected to have equal energy than the initial Ohm mode used
    dB_energy_Ohm=0
    for i in range(rNum-1):
        r=i*dr+dr/2+rmin
        for j in range(thNum-1):
            th=j*dth+dth/2
            dB_energy_Ohm+=(pow(B_field_Ohm[i][j][0]-B_field_k[i][j][0]*sqrt(energy_Ohm/energy),2)
                            +pow(B_field_Ohm[i][j][1]-B_field_k[i][j][1]*sqrt(energy_Ohm/energy),2)
                            +pow(B_field_Ohm[i][j][2]-B_field_k[i][j][2]*sqrt(energy_Ohm/energy),2))*pow(r,2)*sin(th)
    dB_energy_Ohm=dB_energy_Ohm*dth*dr/4

    f.write(str(k*dt) + " " + str(dB_energy) + " " + str(dB_energy_Ohm)+"\n")
    print str(num_file)+" "+str(energy)+" "+str(dB_energy)+" "+str(dB_energy_Ohm)
    k+=plotSteps

f.close()
sys.exit()

