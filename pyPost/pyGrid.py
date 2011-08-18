#!/usr/local/bin/python

# Script:  dump2xyz.py
# Purpose: convert a LAMMPS dump file to XYZ format
# Syntax:  dump2xyz.py dumpfile Nid Ntype Nx Ny Nz xyzfile
#          dumpfile = LAMMPS dump file in native LAMMPS format
#          Nid,Ntype,Nx,Ny,Nz = columns #s for ID,type,x,y,z
#                               (usually 1,2,3,4,5)
#          xyzfile = new XYZ file
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

import sys,os
#path = os.environ["LAMMPS_PYTHON_TOOLS"]
#sys.path.append(path)
from dumpPCOM import dumpPCOM

if len(sys.argv) != 7:
  raise StandardError, "Syntax: pyPost.py dumpfile rootVTK nx ny nz a"

dumpfile = sys.argv[1]
root=sys.argv[2]
#ngrid= [20, 10, 25]
n1 = int (sys.argv[3])
n2 = int (sys.argv[4])
n3 = int (sys.argv[5])
ngrid=[n1, n2, n3]
a=float(sys.argv[6])
d = dumpPCOM(dumpfile)
d.gridProperties(ngrid, root, a)
