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

if len(sys.argv) != 2:
  raise StandardError, "Syntax: pyPost.py dumpfile"

dumpfile = sys.argv[1]
itTypes = ['1.0', '2.0', '3.0']
d = dumpPCOM(dumpfile)
d.writePCOM(itTypes)
