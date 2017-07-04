import re
import sys
import numpy as np
import time

def getCoordinates(inputFile):
    F = open(inputFile,'r')
    step = 1 
    for line in F:
        if "frame" in line:
            s0 = "/home/antadlp/malla/separados-B12N12d-dat/"
            s = s0 + "frame" + str(step) + ".dat"
            line=next(F)                  
            fileOut = open(s, "w")
            for k in range(142): 
                fileOut.write(line.split()[1])
                fileOut.write("   ")
                fileOut.write(line.split()[2])
                fileOut.write("   ")
                fileOut.write(line.split()[3])
                fileOut.write("\n")
                line=next(F)                  

            step += 1 
            fileOut.close()
            
            if step > 5000:
                break

    F.close()
    return;

t = time.time()

inputFile = "./movie-B12N12d.xyz"

getCoordinates(inputFile)

print("elapsed: {}".format(time.time() - t))


