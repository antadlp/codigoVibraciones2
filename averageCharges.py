import pandas as pd
import numpy as np
from get_immediate_subdirectories import *

opt_location = '../optimizations'

dirs = get_immediate_subdirectories(opt_location)

F = open('average-charges.txt', 'w')

for i in dirs:
    xls_location = opt_location + "/" + i + '/charge_mull.xls'
    data = pd.read_csv(xls_location, header=None, delim_whitespace=True)
    atoms = data[1].unique()
    
    print("{}".format(i))
    F.write(i+"\n")
    
    for j in atoms:
        charges = data[data[1].isin([j])].values[:,2]
        print("   Number of {} atoms: {}".format(j,charges.shape[0]))
        F.write("   Number of " + j + " atoms: " + str(charges.shape[0])+"\n")
    
    print("\n")
    F.write("\n")
        
    for j in atoms:
        charges = data[data[1].isin([j])].values[:,2]
        average = charges.sum()/(charges.shape[0])
        print("   Average charge for {}: {}".format(j,average))
        F.write("   Average charge for " + j + ": " + str(average) + "\n")

    
    print("\n")
    F.write("\n")
        
F.close()
