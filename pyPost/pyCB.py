#!/usr/local/bin/python

# Script:  dump2xyz.py
# Purpose: convert a LAMMPS dump file to XYZ format
# Syntax:  dump2xyz.py dumpfile Nid Ntype Nx Ny Nz xyzfile
#          dumpfile = LAMMPS dump file in native LAMMPS format
#          Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
#                               (usually 1,2,3,4,5)
#          xyzfile = new XYZ file
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

import sys,os, glob,  re
#path = os.environ["LAMMPS_PYTHON_TOOLS"]
#sys.path.append(path)

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

    

from dumpPCOM import dumpPCOM

for cual in sys.argv:
    print cual

if len(sys.argv) < 3:
  raise StandardError, "Syntax: pyCB.py input output"
dumpfiles = sys.argv[1]
outputfiles=sys.argv[2]

fdumps=[]
for infile in glob.glob(dumpfiles):   
   fdumps.append(infile)   
   
fs = open(outputfiles, "w")
for fileNa in sorted_nicely(fdumps):       
    if fileNa.endswith('*'):
        continue
    d = dumpPCOM(fileNa)
    d.collisionBIN()
    print>>  fs , d.inicio,  d.final    
fs.close()    
