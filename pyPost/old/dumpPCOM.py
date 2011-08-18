import sys, commands, re, glob, types
from os import popen
from math import *             # any function could be used by set()
import Numeric

#try: from DEFAULTS import PIZZA_GUNZIP
#except: PIZZA_GUNZIP = "gunzip"

# Class definition

class dumpPCOM:

  # --------------------------------------------------------------------

  def __init__(self,*args):
    if len(args) == 0:
        raise StandardError,"dump file name not specified"
    self.snaps = []
    self.nsnaps = self.nselect = 0
    self.names = {}
#    self.tselect = tselect(self)
#    self.aselect = aselect(self)
    self.atype = "type"
    self.bondflag = 0
    self.bondlist = []
    self.triflag = 0
    self.trilist = []
    self.triobj = 0
    self.lineflag = 0
    self.linelist = []

    # flist = list of all dump file names
    
    self.flist = args[0]
    
  def writePCOM(self,*args):    
    if len(args) == 0:
        raise StandardError,"No atom type list specified"    
    idTY = args[0]  
    idx = 1
    nameCOM = []
    fileCOM = []
    nameCOM.append(self.flist+'-AT-'+'all'+'.COM')
    fileCOM.append(open(nameCOM[0],'w'))
    for idtype in idTY:
        nameCOM.append(self.flist+'-AT-'+idtype+'.COM')
        name = nameCOM[idx]
        fpoint = open(nameCOM[idx],'w')
        fileCOM.append(fpoint)
        idx +=1
    self.calcCOM(idTY,fileCOM)
    
  # --------------------------------------------------------------------
  def getCOM(self,snap,idTypes):
      time = snap.time    
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      r = self.names["radius"]
      type = self.names["type"]
      xcom = 0
      ycom = 0
      zcom = 0
      massrT = 0
      atoms = []
      ntypes= len(idTypes)
      xidcom=[]
      yidcom=[]
      zidcom=[]
      idmassrT=[]
      for i in xrange(ntypes):
          xidcom.append(0)
          yidcom.append(0)
          zidcom.append(0)
          idmassrT.append(0)        
      for i in xrange(snap.natoms):            
          atom = snap.atoms[i]
          tipo = atom[type]
          radio = atom[r]
          #massr = atom[r]*atom[r]*atom[r]
          massr = radio*radio*radio
          #atoms.append([atom[id],atom[type],atom[x],atom[y],atom[z]])
          xcom += atom[x]*massr
          ycom += atom[y]*massr
          zcom += atom[z]*massr
          massrT += massr
          for j in xrange(ntypes):
              idT = float(idTypes[j])
              if idT==tipo:
                  xidcom[j]+= atom[x]*massr
                  yidcom[j]+= atom[y]*massr
                  zidcom[j]+= atom[z]*massr
                  idmassrT[j]+= massr
                  break
      if massrT>0:
          xcom/=massrT
          ycom/=massrT
          zcom/=massrT
      respuesta=[]
      respuesta.append(time)
      respuesta.append([xcom, ycom, zcom])          
      for j in xrange(ntypes) :
          if idmassrT[j]>0:
              xidcom[j] /= idmassrT[j]
              yidcom[j] /= idmassrT[j]
              zidcom[j] /= idmassrT[j]
          respuesta.append([xidcom[j], yidcom[j], zidcom[j]])
      return respuesta
  
  def writeVTK(self, *rootF):
      
      if len(*rootF)==0:root="tmp"
      else: root=rootF[0]
      
      file = self.flist
      f=open(file)
      
      snap = self.read_snapshot(f)
      x = self.names["x"]
      y = self.names["y"]
      z = self.names["z"]
      fx= self.names["fx"]
      fy= self.names["fy"]
      fz= self.names["fz"]
      vx= self.names["vx"]
      vy= self.names["vy"]
      vz= self.names["vz"]
      omegax=self.names["omegax"]
      omegay=self.names["omegay"]
      omegaz=self.names["omegaz"]      
      r = self.names["radius"]      
      type = self.names["type"]
      try:
        epot=self.names["c_epot_n"]
        epot_flag=1
      except:
         epot_flag=0 
      try:
        ekin=self.names["c_enKin"]
        ekin_flag=1
      except:
         ekin_flag=0          
      try:
        cn=self.names["f_cn"]
        cn_flag=1
      except:
         cn_flag=0                
      n=0
      while snap:
          vtkFile = root+str(n).zfill(10)+".vtk"
          vtkFile_bb = root+str(n).zfill(10)+"_boundingBox.vtk"
          print vtkFile, vtkFile_bb
          time = snap.time
          xlo=snap.xlo
          xhi=snap.xhi
          ylo=snap.ylo
          yhi=snap.yhi
          zlo=snap.zlo
          zhi=snap.zhi
          self.boundingBox(vtkFile_bb,xlo,xhi,ylo,yhi,zlo,zhi)
          
          atoms=snap.atoms
          
          fs = open(vtkFile,"w")  
          print >>fs,"# vtk DataFile Version 2.0"
          print >>fs,"Generated by pizza.py"
          print >>fs,"ASCII"
          print >>fs,"DATASET POLYDATA"
          print >>fs,"POINTS %d float" % len(atoms)
          for atom in atoms:
              print >>fs,atom[x],atom[y],atom[z]  #write x,y,z  [atom[0]=id, atom[1]=type]
          print >>fs,"VERTICES",len(atoms),2*len(atoms)
          for i in xrange(len(atoms)):
              print >>fs,1,i
          print >>fs,"POINT_DATA",len(atoms)

          print >>fs,"VECTORS ","f"," float"
          for atom in atoms:
              print >>fs,atom[fx],atom[fy],atom[fz]
          print >>fs,"VECTORS ","v"," float"
          for atom in atoms:
              print >>fs,atom[vx],atom[vy],atom[vz]                        
          print >>fs,"VECTORS ","omega"," float"
          for atom in atoms:
              print >>fs,atom[omegax],atom[omegay],atom[omegaz]              
     
          print >>fs,"SCALARS atom_type int 1"
          print >>fs,"LOOKUP_TABLE default"
          for atom in atoms:  #loop all atoms
              itype = int(atom[type])
              print >>fs,itype,

          print >>fs,"SCALARS radius float 1"
          print >>fs,"LOOKUP_TABLE default"
          for atom in atoms:  #loop all atoms
              print >>fs,atom[r]  
          if epot_flag:
              print >>fs,"SCALARS epot float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[epot]  
          if ekin_flag:
              print >>fs,"SCALARS ekin float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[ekin]                  
          if cn_flag:
              print >>fs,"SCALARS CoordN float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cn]                                  
          fs.close()          
          snap = self.read_snapshot(f)
          n+=1
          
          
          
  def boundingBox(self,file,xlo,xhi,ylo,yhi,zlo,zhi):
      f = open(file,"w")  
      print >>f,"# vtk DataFile Version 2.0"
      print >>f,"Generated by pizza.py"
      print >>f,"ASCII"
      print >>f,"DATASET RECTILINEAR_GRID"
      print >>f,"DIMENSIONS 2 2 2"  
      print >>f,"X_COORDINATES 2 float"
      print >>f,xlo,xhi
      print >>f,"Y_COORDINATES 2 float"
      print >>f,ylo,yhi  
      print >>f,"Z_COORDINATES 2 float"
      print >>f,zlo,zhi
 
 
  def calcCOM(self, idTypes, filePointer):
  
    file=self.flist
    f=open(file)

#    snap = self.read_snapshot(f)
#    comCoords = self.getCOM(snap, idTypes)
#    time = comCoords[0]
#    for i in xrange(len(filePointer)):
#        cual = comCoords[i+1]
#        print>>filePointer[i],time,cual[0], cual[1], cual[2]
    
    snap = self.read_snapshot(f)
    while snap:
        print snap.time,
        sys.stdout.flush()
        comCoords = self.getCOM(snap, idTypes)
        time = comCoords[0]
        for i in xrange(len(filePointer)):
            cual = comCoords[i+1]
            print>>filePointer[i],time,cual[0], cual[1], cual[2]
        snap = self.read_snapshot(f)
                    
    for i in xrange(len(filePointer)):
        filePointer[i].close()
    print


  
  def read_snapshot(self,f):
    try:
      snap = Snap()
      item = f.readline()
      snap.time = int(f.readline().split()[0])    # just grab 1st field
      item = f.readline()
      snap.natoms = int(f.readline())

      snap.aselect = Numeric.zeros(snap.natoms)

      item = f.readline()
      words = f.readline().split()
      snap.xlo,snap.xhi = float(words[0]),float(words[1])
      words = f.readline().split()
      snap.ylo,snap.yhi = float(words[0]),float(words[1])
      words = f.readline().split()
      snap.zlo,snap.zhi = float(words[0]),float(words[1])

      item = f.readline()
      if len(self.names) == 0:
        words = item.split()[2:]
        if len(words):
          for i in range(len(words)):
            if words[i] == "xs" or words[i] == "xu":
              self.names["x"] = i
            elif words[i] == "ys" or words[i] == "yu":
              self.names["y"] = i
            elif words[i] == "zs" or words[i] == "zu":
              self.names["z"] = i
            else: self.names[words[i]] = i

      if snap.natoms:
        words = f.readline().split()
        ncol = len(words)
        for i in xrange(1,snap.natoms):
          words += f.readline().split()
        floats = map(float,words)
        atoms = Numeric.zeros((snap.natoms,ncol),Numeric.Float)
        start = 0
        stop = ncol
        for i in xrange(snap.natoms):
          atoms[i] = floats[start:stop]
          start = stop
          stop += ncol
      else: atoms = None
      snap.atoms = atoms
      return snap
    except:
      return 0

## --------------------------------------------------------------------
## one snapshot

class Snap:
  pass
