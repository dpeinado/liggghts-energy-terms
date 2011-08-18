import sys, commands, re, glob, types
from os import popen
from math import *             # any function could be used by set()
import Numeric
from scitools.std import *

class dumpThermo:  

  def __init__(self,*args):
    if len(args) == 0:
        raise StandardError,"dump file name not specified"
    self.snaps = []
    self.names = []
    self.flist = args[0]
  def plotOutput(self,  *rootF):
    if len(*rootF)==0:root="tmp"
    else: root=rootF[0]
    file = root+"01.dat"
    f=open(file)
    names=f.readline().split()
    ncol=len(names)
    flag = 1
    allwords = []
    ncount=0
    while flag:
        line=f.readline()
        if line=="":
            break
        words = line.split()
        allwords += words        
        ncount+=1
    floats = map(float,allwords)   
    atoms = Numeric.zeros((ncount,ncol),Numeric.Float)
    start = 0
    stop = ncol
    for i in xrange(ncount):
        atoms[i] = floats[start:stop]
        start = stop
        stop += ncol
    x = atoms[:, 0]
    y = atoms[:, 2]
   #plot(atoms[:, 0], atoms[:, 2], 'y-')
   #hold('on')
   #plot(atoms[:, 0], atoms[:, 3], 'm-')
   #plot(atoms[:, 0], atoms[:, 4], 'c-')
   #plot(atoms[:, 0], atoms[:, 5], 'r-')
   #plot(atoms[:, 0], atoms[:, 6], 'g-')
   #plot(atoms[:, 0], atoms[:, 7], 'b-')
   #plot(atoms[:, 0], atoms[:, 8], 'k-')
   #plot(atoms[:, 0], atoms[:, 9], 'b--')
   #plot(atoms[:, 0], atoms[:, 10], 'r--')
    Energy0=Numeric.zeros( (len(atoms[:, 1]), 1), Numeric.Float)    
    E0 = atoms[0, 2]
    Energy0[0:len(Energy0)] = E0
    plot(atoms[:, 0], atoms[:, 2], 'y-',legend=names[2],linewidth=2.0)
    hold('on')
    plot(atoms[:, 0], atoms[:, 3], 'm-',legend=names[3],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 4], 'c-',legend=names[4],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 5], 'r-',legend=names[5],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 6], 'g-',legend=names[6],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 7], 'b-',legend=names[7],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 8], 'k-',legend=names[8],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 9], 'b-',legend=names[9],linewidth=2.0)
    plot(atoms[:, 0], atoms[:, 10], 'r-',legend=names[10],linewidth=2.0)
    plot(atoms[:, 0],  Energy0, 'r-', legend='E0')
    grid(True)
    
    #legend(names[2],  names[3],  names[4], names[5],  names[6],  names[7], names[8],  names[9], names[10])
    raw_input("hola")
    
  def thermoOutput(self, *rootF):
      
      if len(*rootF)==0:root="tmp"
      else: root=rootF[0]
      
      file = self.flist
      f=open(file)
      flag=1
      n=0
      while flag:
        snap = self.read_snapshot(f)
        if not snap:
            break
        vtkFile = root+str(n).zfill(2)+".dat"
        fs = open(vtkFile,"w")
        atoms = snap.atoms
        ncol = len(self.names)
        linea = ""
        for name in self.names:
            linea += name + " "
        print >> fs,  linea
        for atom in atoms:
            linea = ""
            for i in xrange(ncol):
                linea += repr(atom[i]) + " "
            print >> fs,  linea
        fs.close()
        n+=1
  def read_snapshot(self,f):
    try:
        snap = Snap()
        flag = 1
        while flag:
            line = f.readline()
            if line =="":
                return 0
            words = line.split()            
            if len(words)>0:
                if (words[0] == "Setting") and (words[2] == "run"):
                    flag = 0
        line = f.readline()
        self.names = f.readline().split()
        ncol = len(self.names)
        flag = 1
        allwords = []
        snap.nlines = 0
        while flag:
            words = f.readline().split()
            if words[0] == "Loop":
                flag = 0
            else:
                allwords += words
                snap.nlines+=1
        floats = map(float,allwords)   
        atoms = Numeric.zeros((snap.nlines,ncol),Numeric.Float)
        start = 0
        stop = ncol
        for i in xrange(snap.nlines):
          atoms[i] = floats[start:stop]
          start = stop
          stop += ncol        
        snap.atoms = atoms
        return snap
    except:
        return 0      
  def column(matrix, i):
    return [row[i] for row in matrix]
            
class Snap:
  pass
          
