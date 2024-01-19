#
#  Common Python functions used in multiple notebooks
#
from thermoengine import core
import numpy as np
#
def print_species(n, s):
    fmt = '{0:>10s} {1:13.6g}'
    print (fmt.format('      SiO2', n[0]+s[0]+s[1]+3*s[2]))
    print (fmt.format('      TiO2', n[1]))
    print (fmt.format('     Al2O3', n[2]+2*s[2]))
    print (fmt.format('     Fe2O3', n[3]))
    print (fmt.format('     Cr2O3', n[4]))
    print (fmt.format('   Fe2SiO4', n[5]))
    print (fmt.format(' MnSi1/2O2', n[6]))
    print (fmt.format('   Mg2SiO4', n[7]))
    print (fmt.format(' NiSi1/2O2', n[8]))
    print (fmt.format(' CoSi1/2O2', n[9]))
    print (fmt.format('    CaSiO3', n[10]-s[0]))
    print (fmt.format('   Na2SiO3', n[11]-2*s[1]))
    print (fmt.format('   KAlSiO4', n[12]-4*s[2]))
    print (fmt.format(' Ca3(PO4)2', n[13]))
    print (fmt.format('       H2O', n[14]))
    print (fmt.format('       CO2', n[15]-s[0]))
    print (fmt.format('    ZrSiO4', n[16]-s[1]-s[2]))
    print (fmt.format('     CaCO3', s[0]))
    print (fmt.format('Na4ZrSi2O8', s[1]))
    print (fmt.format(' K4ZrSi2O8', s[2]))
    return
