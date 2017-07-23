import pandas as pd
import numpy as np 
import os, sys
import time



def getXYZFromMovie(inputFile,steps,atoms):

    F = open(inputFile,'r')

    xList =[]
    yList =[]
    zList =[]
    step = 1 
    flag1 = True

    while flag1 == True:
        line=next(F)
        if "frame" in line:
            k=1 #is the counter for the index of every atom
            j=1 #counter for a constant empty lines
            x = []
            y = []
            z = []
            Za = []
            for j in range(1):#empty constant lines
                line = next(F)
            
            l=1
            for k in range(atoms): 
                Za.append(line.split()[0])
                x.append(float(line.split()[1])) #get x coordinate
                y.append(float(line.split()[2])) #get y coordinate
                z.append(float(line.split()[3])) #get z coordinate
                if ( ((step + 1) > steps) & (l == atoms) ):
                    flag1=False
                else:
                    l+=1
                    line=next(F) #jump to next atom

            step += 1 #when this loops ends, the outside loop search for the
                      #next step on the dynamic that's why step+=1
            xList.append(x)
            yList.append(y)
            zList.append(z)


    nx = np.array(xList)
    ny = np.array(yList)
    nz = np.array(zList)
    F.close()
    
    return nx, ny, nz, Za;

