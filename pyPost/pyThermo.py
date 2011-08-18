import sys,os
import subprocess
from dumpThermo import dumpThermo

#if len(sys.argv) < 2:
 # raise StandardError, "Syntax: pyPost.py dumpfile [vtkROOT]"

#dumpfile = sys.argv[1]
dumpfile = "output.log"
#if len(sys.argv)==3:
#   vtkROOT = sys.argv[2]
#else: vtkROOT = "tmp"
d = dumpThermo(dumpfile)
d.thermoOutput("output")
d.plotOutput("output")
