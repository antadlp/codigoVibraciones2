import pandas as pd
import numpy as np 
import os, sys
import time



def getXYZFromXYZF(inputFile, atoms):

    F = open(inputFile,'r')

    x = []
    y = []
    z = []
    Za = []

    line = next(F)
    line = next(F)
            
    for line in F: 
        Za.append(line.split()[0])
        x.append(float(line.split()[1])) #get x coordinate
        y.append(float(line.split()[2])) #get y coordinate
        z.append(float(line.split()[3])) #get z coordinate

    nx = np.array(x)
    ny = np.array(y)
    nz = np.array(z)

    F.close()
    
    return nx, ny, nz, Za;

