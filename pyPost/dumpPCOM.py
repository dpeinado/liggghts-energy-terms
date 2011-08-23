import sys, commands, re, glob, types
from os import popen
from math import *             # any function could be used by set()
#import Numeric
from numpy import *
#from evtk.hl import gridToVTK 

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
#        self.tselect = tselect(self)
#       self.aselect = aselect(self)
        self.atype = "type"
        self.bondflag = 0
        self.bondlist = []
        self.triflag = 0
        self.trilist = []
        self.triobj = 0
        self.lineflag = 0
        self.linelist = []
#       flist = list of all dump file names    
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
        cpen=self.names["f_CPEn"]
        cpen_flag=1
      except:
         cpen_flag=0 
      try:
        cden=self.names["f_CDEn"]
        cden_flag=1
      except:
         cden_flag=0          
      try:
        cpet=self.names["f_CPEt"]
        cpet_flag=1
      except:
         cpet_flag=0             
      try:
        cdetv=self.names["f_CDEVt"]
        cdetv_flag=1
      except:
         cdetv_flag=0          
      try:
        cdetf=self.names["f_CDEFt"]
        cdetf_flag=1
      except:
         cdetf_flag=0                   
      try:
        ctfw=self.names["f_CTFW"]
        ctfw_flag=1
      except:
         ctfw_flag=0                  
      try:
        deh=self.names["f_DEH"]
        deh_flag=1
      except:
         deh_flag=0                       
      try:
          mass=self.names["mass"]
          mass_flag=1
      except:
          mass_flag=0          
      try:
          ePGp=self.names["v_ePGp"]
          ePGp_flag=1
      except:
          ePGp_flag=0       
      try:
          eKinLp=self.names["c_eKinLp"]
          eKinLp_flag=1
      except:
          eKinLp_flag=0              
      try:
          eKinRp=self.names["c_eKinRp"]
          eKinRp_flag=1
      except:
          eKinRp_flag=0   
                             
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
          print 'Natoms = ', len(atoms)
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
          if cpen_flag:
              print >>fs,"SCALARS CPEn float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cpen]  
          if cden_flag:
              print >>fs,"SCALARS CDEn float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cden]  
          if cpet_flag:
              print >>fs,"SCALARS CPEt float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cpet]  
          if cdetv_flag:
              print >>fs,"SCALARS CDEVt float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cdetv]  
          if cdetf_flag:
              print >>fs,"SCALARS CDEFt float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[cdetf]                  
          if ctfw_flag:
              print >>fs,"SCALARS CTFW float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[ctfw]  
          if deh_flag:
              print >>fs,"SCALARS DEH float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[deh]                              
          if mass_flag:
              print >>fs,"SCALARS Mass float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[mass]                                              
          if ePGp_flag:
              print >>fs,"SCALARS ePGp float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[ePGp]        
          if eKinLp_flag:
              print >>fs,"SCALARS eKinLp float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[eKinLp]                      
          if eKinRp_flag:
              print >>fs,"SCALARS eKinRp float 1"
              print >>fs,"LOOKUP_TABLE default"
              for atom in atoms:  #loop all atoms
                print >>fs,atom[eKinRp]                                                      
          fs.close()          
          snap = self.read_snapshot(f)
          n+=1
    def collisionBIN(self):
      firstTime=True
      file = self.flist
      f=open(file)      
      snap = self.read_snapshot(f)
      n=snap.natoms
      if n <> 2:
        raise StandardError,"Number of atoms <> 2"    
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
         
      while snap:
        atoms=snap.atoms
        atom1=atoms[0]
        atom2=atoms[1]
        delx = atom1[x] - atom2[x]
        dely = atom1[y] - atom2[y]
        delz = atom1[z] - atom2[z]
        rsq = delx*delx+dely*dely+delz*delz
        rmod=sqrt(rsq)
        r1 = atom1[r]
        r2 = atom2[r]
        radsum=r1+r2
        deltan=radsum-rmod
      
        vr1 = atom1[vx]-atom2[vx]
        vr2 = atom1[vy]-atom2[vy]
        vr3 = atom1[vz]-atom2[vz]
        vnnr = (vr1*delx+vr2*dely+vr3*delz)/rmod
      
        vn1 = delx*vnnr/rsq
        vn2 = dely*vnnr/rsq
        vn3 = delz*vnnr/rsq
      
        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;
        cr1 = r1-0.5*deltan;
        cr2 = r2-0.5*deltan;
        wr1 = (cr1*atom1[omegax] + cr2*atom2[omegax]) /rmod;
        wr2 = (cr1*atom1[omegay] + cr2*atom2[omegay]) /rmod;
        wr3 = (cr1*atom1[omegaz] + cr2*atom2[omegaz]) /rmod;
        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);                
        if firstTime:
            phi = fabs(vrel/vnnr)
            self.inicio=phi
            self.vtr1=vtr1
            self.vtr2=vtr2
            self.vtr3=vtr3
            self.vni=vnnr
            firstTime=False
        else:
            phi = fabs(vrel/self.vni)
        modulo = self.vtr1*vtr1+self.vtr2*vtr2+self.vtr3*vtr3
        if modulo <>0.0:
            signo=modulo/fabs(modulo)            
        else:
            signo=0.0
        phi=phi*signo
        snap = self.read_snapshot(f)
      self.final = phi
    def writeGRAPH(self, root, enFile):      
        energyFlag = False
        IKE = 0.0
        if enFile <> None:
            energyFlag= True
            fen = open(enFile)        
        file = self.flist
        f=open(file)
      
        snap = self.read_snapshot(f)
        n=snap.natoms
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
        cpen=self.names["f_CPEn"]
        cden=self.names["f_CDEn"]
        cpet=self.names["f_CPEt"]
        cdevt=self.names["f_CDEVt"]
        cdeft=self.names["f_CDEFt"]
        ctfw=self.names["f_CTFW"]
        deh=self.names["f_DEH"]
        mass=self.names["mass"]

        atomFile=[]
        fs=[]
        strtmp=""
        for i in xrange(n+1):
            atomFile.append(root+str(i)+'.dat')
            print atomFile[i]
            fs.append(open(atomFile[i], "w"))
        if energyFlag:
            line=  fen.readline()
            line=  fen.readline()        
    
        while snap:
            if energyFlag:
                line=  fen.readline()
                fields = line.split()
                IKE = float(fields[1])
            time = snap.time                  
            atoms=snap.atoms
            atom=atoms[0]
            kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
            kR = 0.5*2.0/5.0*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
            kET = kE
            kRT = kR
            cpenT = atom[cpen]
            cpetT = atom[cpet]
            cdenT = atom[cden]
            cdetTV = atom[cdevt]
            cdetTF = atom[cdeft]
            ctfwT = atom[ctfw]
            dehT = atom[deh]

            EnCons = kE+kR+atom[cpen]
            EnTot = EnCons+atom[ctfw]+atom[deh]+atom[cden]+atom[cdevt]+atom[cdeft]
            EnConsT = EnCons
            EnTotT = EnTot
            print >>fs[0], time, kE,  kR,  atom[cpen], atom[cpet],  atom[cden], atom[cdevt], atom[cdeft], atom[ctfw], atom[deh], EnCons, EnTot
            for i in xrange(1, n):
                atom=atoms[i]
                kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
                kR = 0.5*2.0/5.0*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
                kET += kE
                kRT += kR
                cpenT += atom[cpen]
                cpetT += atom[cpet]
                cdenT += atom[cden]
                cdetTV += atom[cdevt]
                cdetTF += atom[cdeft]
                ctfwT += atom[ctfw]
                dehT += atom[deh]          

                EnCons = kE+kR+atom[cpen]
                EnTot = EnCons+atom[ctfw]+atom[deh]+atom[cden]+atom[cdevt]+atom[cdeft]
                EnConsT += EnCons
                EnTotT += EnTot
                print >>fs[i], time, kE,  kR,  atom[cpen], atom[cpet],  atom[cden], atom[cdevt], atom[cdeft], atom[ctfw], atom[deh], EnCons, EnTot
            print >>fs[n],  time, kET, kRT, cpenT, cpetT, cdenT, cdetTV, cdetTF, ctfwT, dehT, EnConsT, (EnTotT-IKE)
            snap = self.read_snapshot(f)
        for i in xrange(n+1):
            fs[i].close()
        
    def checkEnergy(self, root, enFile):      
        energyFlag = False
        IKE = 0.0
        if enFile <> None:
            print "Energy Flag set ON"
            energyFlag= True
            fen = open(enFile)        
        file = self.flist
        f=open(file)
        
        snap = self.read_snapshot(f)
        n=snap.natoms
        if n>1:
            raise StandardError, "Solo quiero una particula (energia particula contra pared)"
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
        cpen=self.names["f_CPEn"]
        cden=self.names["f_CDEn"]
        cpet=self.names["f_CPEt"]
        cdevt=self.names["f_CDEVt"]
        cdeft=self.names["f_CDEFt"]
        ctfw=self.names["f_CTFW"]
        deh=self.names["f_DEH"]
        mass=self.names["mass"]
        try:
            epg=self.names["v_ePGp"]
            epg_flag=1
            print "Gravitational Potential Energy Flag set ON"
        except:
            epg_flag=0 
        try:
            ekinl=self.names["c_eKinLp"]
            ekinl_flag=1
            print "Translational Kinetic Energy Flag set ON"
        except:
            ekinl_flag=0          
        try:
            ekinr=self.names["v_eKinRp"]
            ekinr_flag=1
            print "Rotational Kinetic Energy Flag set ON"
        except:
            ekinr_flag=0          
            
        atomFile=root+'.dat'    
        strtmp=""
        fs = open(atomFile, "w")
        if energyFlag:
            line=  fen.readline()
            line=  fen.readline()        
        print >>fs, '# time', 'eKIN',  'eCOL',  'ePG',  'DEH',  'IKE', 'LHS',  'ErrorA',  'ErrorR',   'kE', 'kR',  'eKIN',  'cpeN',  'cpeT',  'cdeN',  'cdeTV',  'cdeTF',  'ctfW'
        n=0
        while snap:
            if energyFlag:
                line=  fen.readline()
                fields = line.split()
                IKE = float(fields[1])
            time = snap.time                  
            atoms=snap.atoms
            atom=atoms[0]
            if ekinl_flag:
                kE = atom[ekinl]
            else:  kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
            if ekinr_flag:
                kR = atom[ekinr]
            else: kR = 0.2*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
            eKIN = kE+kR
            cpeN = atom[cpen]
            cpeT = atom[cpet]
            cdeN = atom[cden]
            cdeTV = atom[cdevt]
            cdeTF = atom[cdeft]
            ctfW = atom[ctfw]
            eCOL= cpeN+cpeT+cdeN+cdeTV+cdeTF+ctfW
            DEH = atom[deh]
            if(epg_flag):
                ePG = atom[epg]
            else: ePG = 0.0
            if(n==0):
                En0=eKIN+eCOL+ePG
            LHS = eKIN+eCOL+ePG+DEH
            RHS = IKE+En0
            ErrorA = math.fabs(LHS-RHS)
            Denom=0.5*(LHS+RHS)
            if fabs(Denom)>0.0:
                ErrorR = ErrorA/Denom*100.0
            else: ErrorR =0
            print >>fs, time, eKIN,  eCOL,  ePG,  DEH,  IKE,  LHS,  ErrorA,  ErrorR, kE,  kR,  eKIN,  cpeN,  cpeT,  cdeN,  cdeTV,  cdeTF,  ctfW
            snap = self.read_snapshot(f)        
            n+=1
        fs.close()
            
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

            snap.aselect = zeros(snap.natoms)

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
                atoms = zeros((snap.natoms,ncol), dtype=float32)
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
    def gridProperties(self, *args):
        if len(args) == 0:
            raise StandardError,"No atom type list specified"    
        nx=args[0][0]
        ny=args[0][1]
        nz=args[0][2]
        npx=nx+1
        npy=ny+1
        npz=nz+1
        root=args[1]
        na=args[2]
        file=self.flist
        f=open(file)
        snap = self.read_snapshot(f)
        n=snap.natoms
        xlo=snap.xlo
        xhi=snap.xhi
        ylo=snap.ylo
        yhi=snap.yhi
        zlo=snap.zlo
        zhi=snap.zhi
        lx=xhi-xlo
        ly=yhi-ylo
        lz=zhi-zlo
        dx=lx/nx
        dy=ly/ny
        dz=lz/nz
        print 'Simulation box: ',  xlo,  xhi,  ylo,  yhi,  zlo,  zhi
        print 'Grid discretization: ', nx, ny, nz
        print 'Grid spacing: ',  dx, dy, dz
        print 'Cube dimension %.2f times diameter'% na
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
        cpen=self.names["f_CPEn"]
        cden=self.names["f_CDEn"]
        cpet=self.names["f_CPEt"]
        cdevt=self.names["f_CDEVt"]
        cdeft=self.names["f_CDEFt"]
        ctfw=self.names["f_CTFW"]
        deh=self.names["f_DEH"]
        mass=self.names["mass"]
        
        nfile=0        
        while snap:
            #massTotal1=0
            #massTotal2=0
            cpenEuler=zeros( (nx, ny, nz, 4) )
            cdenEuler=zeros( (nx, ny, nz, 4) )
            cdevtEuler=zeros( (nx, ny, nz, 4) )
            cdeftEuler=zeros( (nx, ny, nz, 4) )
            dehEuler=zeros( (nx, ny, nz, 4) )
            massEuler=zeros( (nx, ny, nz, 4) )
            epgEuler = zeros( (nx, ny, nz, 4) )
            ekEuler = zeros( (nx, ny, nz, 4) )
            
            atoms=snap.atoms
            for atom in atoms:
                a=2*na*atom[r]
                #print "na = %s radius = %s a = %s" %(na, atom[r], a)
                #massTotal1+=atom[mass]
# calculation of volume factor for x dimension
                ncx=int(a/dx)+2
                if ncx%2==0:
                    ncx+=1
                xindexes=zeros(ncx)
                fxcell=zeros(ncx)
                indlo=-int(ncx/2)
                i=int((atom[x]-xlo)/dx)
                xcubelo=atom[x]-a/2
                if(xcubelo<xlo):
                    xcubelo+=lx
                xcubehi=atom[x]+a/2
                if(xcubehi>xhi):
                    xcubehi-=lx
                for m in range(ncx):
                    xindexes[m]=i+(indlo+m)
                    if(xindexes[m]<0):
                        xindexes[m]+=nx
                    if(xindexes[m]>=nx):
                        xindexes[m]-=nx
                for m in range(ncx):
                    xcelllo=xindexes[m]*dx+xlo
                    xcellhi=(xindexes[m]+1)*dx+xlo
                    if((xcubelo>xcelllo)&(xcubelo<xcellhi)):
                        indexCubeLo=m
                    if((xcubehi>xcelllo)&(xcubehi<xcellhi)):
                        indexCubeHi=m
                for m in range(indexCubeLo):
                    fxcell[m]=0
                m=indexCubeLo
                xcelllo=xindexes[m]*dx+xlo
                xcellhi=(xindexes[m]+1)*dx+xlo                   
                fxcell[indexCubeLo]=(xcellhi-xcubelo)/dx
                for m in range(indexCubeLo+1, indexCubeHi):
                    fxcell[m]=1
                m=indexCubeHi
                xcelllo=xindexes[m]*dx+xlo
                xcellhi=(xindexes[m]+1)*dx+xlo                   
                fxcell[indexCubeHi]=(xcubehi-xcelllo)/dx
                for m in range(indexCubeHi+1, ncx):
                    fxcell[m]=0
# calculation of volume factor for y dimension
                ncy=int(a/dy)+2
                if ncy%2==0:
                    ncy+=1
                yindexes=zeros(ncy)
                fycell=zeros(ncy)
                indlo=-int(ncy/2)
                j=int((atom[y]-ylo)/dy)
                ycubelo=atom[y]-a/2
                if(ycubelo<ylo):
                    ycubelo+=ly
                ycubehi=atom[y]+a/2
                if(ycubehi>yhi):
                    ycubehi-=ly
                for m in range(ncy):
                    yindexes[m]=j+(indlo+m)
                    if(yindexes[m]<0):
                        yindexes[m]+=ny
                    if(yindexes[m]>=ny):
                        yindexes[m]-=ny
                for m in range(ncy):
                    ycelllo=yindexes[m]*dy+ylo
                    ycellhi=(yindexes[m]+1)*dy+ylo
                    if((ycubelo>ycelllo)&(ycubelo<ycellhi)):
                        indexCubeLo=m
                    if((ycubehi>ycelllo)&(ycubehi<ycellhi)):
                        indexCubeHi=m
                for m in range(indexCubeLo):
                    fycell[m]=0
                m=indexCubeLo
                ycelllo=yindexes[m]*dy+ylo
                ycellhi=(yindexes[m]+1)*dy+ylo
                fycell[indexCubeLo]=(ycellhi-ycubelo)/dy
                for m in range(indexCubeLo+1, indexCubeHi):
                    fycell[m]=1
                m=indexCubeHi
                ycelllo=yindexes[m]*dy+ylo
                ycellhi=(yindexes[m]+1)*dy+ylo
                fycell[indexCubeHi]=(ycubehi-ycelllo)/dy
                for m in range(indexCubeHi+1, ncy):
                    fycell[m]=0                    
 # calculation of volume factor for z dimension
                ncz=int(a/dz)+2
                if ncz%2==0:
                    ncz+=1
                zindexes=zeros(ncz)
                fzcell=zeros(ncz)
                indlo=-int(ncz/2)
                k=int((atom[z]-zlo)/dz)
                zcubelo=atom[z]-a/2
                if(zcubelo<zlo):
                    zcubelo+=lz
                zcubehi=atom[z]+a/2
                if(zcubehi>zhi):
                    zcubehi-=lz
                for m in range(ncz):
                    zindexes[m]=k+(indlo+m)
                    if(zindexes[m]<0):
                        zindexes[m]+=nz
                    if(zindexes[m]>=nz):
                        zindexes[m]-=nz
                for m in range(ncz):
                    zcelllo=zindexes[m]*dz+zlo
                    zcellhi=(zindexes[m]+1)*dz+zlo
                    if((zcubelo>zcelllo)&(zcubelo<zcellhi)):
                        indexCubeLo=m
                    if((zcubehi>zcelllo)&(zcubehi<zcellhi)):
                        indexCubeHi=m
                for m in range(indexCubeLo):
                    fzcell[m]=0
                m=indexCubeLo
                zcelllo=zindexes[m]*dz+zlo
                zcellhi=(zindexes[m]+1)*dz+zlo
                fzcell[indexCubeLo]=(zcellhi-zcubelo)/dz
                for m in range(indexCubeLo+1, indexCubeHi):
                    fzcell[m]=1
                m=indexCubeHi
                zcelllo=zindexes[m]*dz+zlo
                zcellhi=(zindexes[m]+1)*dz+zlo
                fzcell[indexCubeHi]=(zcubehi-zcelllo)/dz
                for m in range(indexCubeHi+1, ncz):
                    fzcell[m]=0
                Vcell=dx*dy*dz
                Vcube=a*a*a
                #fraccionTotal=0
                #print "a = %s Vcell = %s Vcube = %s" %(a, Vcell, Vcube)
                #massEuler=zeros( (nx, ny, nz) )
                for i in range(ncx):
                    for j in range(ncy):
                        for k in range(ncz):
                            fractionVol=fxcell[i]*fycell[j]*fzcell[k]*Vcell/Vcube
                            #print "fractionVol = %s" %fractionVol
                            #fraccionTotal+=fractionVol
                            massEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*atom[mass]
                            cpenEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*atom[cpen]
                            dehEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]]   +=fractionVol*atom[deh]
                            cdenEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*atom[cden]
                            cdevtEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]]+=fractionVol*atom[cdevt]
                            cdeftEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*atom[cdeft]                              
                            epgEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*atom[mass]*atom[z]
                            kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
                            kR = 0.5*2.0/5.0*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
                            ekEuler[xindexes[i], yindexes[j], zindexes[k], atom[type]] +=fractionVol*(kE+kR)
                            
#                print "fraccion Total = %s" %fraccionTotal
#                masaporahora=0
#                for i in range(nx):
#                    for j in range(ny):
#                        for k in range(nz):
#                            masaporahora+=massEuler[i, j, k]
#                print "Masa por ahora = %s Masa particula = %s" %(masaporahora, atom[mass])
            vtkFile = root+str(nfile).zfill(10)+".vtk"
            fs = open(vtkFile,"w")  
            DissipationEuler=zeros( (nx, ny, nz, 4) )
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        massEuler[i, j, k, 0]   = massEuler[i, j, k, 1] + massEuler[i, j, k, 2]+ massEuler[i, j, k, 3]
                        cpenEuler[i, j, k, 0]   = cpenEuler[i, j, k, 1] + cpenEuler[i, j, k, 2]+ cpenEuler[i, j, k, 3]
                        dehEuler[i, j, k, 0]     = dehEuler[i, j, k, 1] + dehEuler[i, j, k, 2]+ dehEuler[i, j, k, 3]
                        cdenEuler[i, j, k, 0]   = cdenEuler[i, j, k, 1] + cdenEuler[i, j, k, 2]+ cdenEuler[i, j, k, 3]
                        cdevtEuler[i, j, k, 0]  = cdevtEuler[i, j, k, 1]+ cdevtEuler[i, j, k, 2]+ cdevtEuler[i, j, k, 3]
                        cdeftEuler[i, j, k, 0]   = cdeftEuler[i, j, k, 1] + cdeftEuler[i, j, k, 2]+ cdeftEuler[i, j, k, 3]
                        epgEuler[i, j, k, 0]     = epgEuler[i, j, k, 1]+ epgEuler[i, j, k, 2]+ epgEuler[i, j, k, 3]
                        ekEuler[i, j, k, 0]       = ekEuler[i, j, k, 1]+ ekEuler[i, j, k, 2]+ ekEuler[i, j, k, 3]
                        for l in range(4):
                            DissipationEuler[i, j, k, l]=dehEuler[i, j, k, l]+cdenEuler[i, j, k, l]+cdeftEuler[i, j, k, l]+cdevtEuler[i, j, k, l]
                        #massTotal2+=massEuler[i, j, k]
            #print "Masa particulas = %s Masa Euler =%s" % (massTotal1,  massTotal2)
            print "Writing %s" % vtkFile
            print >>fs,"# vtk DataFile Version 3.0"
            print >>fs,"Generated by pizza.py"
            print >>fs,"ASCII"
            print >>fs,"DATASET RECTILINEAR_GRID"
            print >>fs,"DIMENSIONS" , npx, npy, npz
            print >>fs, "X_COORDINATES" ,  npx,  "double"
            for ii in range(npx):
                print >>fs, xlo+ii*dx 
            print >>fs, "Y_COORDINATES" ,  npy,  "double"
            for ii in range(npy):
                print >>fs, ylo+ii*dy                 
            print >>fs, "Z_COORDINATES" ,  npz,  "double"
            for ii in range(npz):
                print >>fs, zlo+ii*dz
            totalNumber=nx*ny*nz
            print >>fs,"CELL_DATA",totalNumber
            print >>fs, "FIELD DEHdata 1"            
            print >>fs,'DEH 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,dehEuler[i, j, k, 0], dehEuler[i, j, k, 1], dehEuler[i, j, k, 2], dehEuler[i, j, k, 3]
            print >>fs, "FIELD CDENdata 1"                                    
            print >>fs,'CDEn 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,cdenEuler[i, j, k, 0], cdenEuler[i, j, k, 1], cdenEuler[i, j, k, 2], cdenEuler[i, j, k, 3]                    
            print >>fs, "FIELD CDEVTdata 1"                                    
            print >>fs,'CDEVt 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs, cdevtEuler[i, j, k, 0], cdevtEuler[i, j, k, 1], cdevtEuler[i, j, k, 2], cdevtEuler[i, j, k, 3]               
            print >>fs, "FIELD CDEFTdata 1"                                    
            print >>fs,'CDEFt 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,cdeftEuler[i, j, k, 0], cdeftEuler[i, j, k, 1], cdeftEuler[i, j, k, 2], cdeftEuler[i, j, k, 3]                                        
            print >>fs, "FIELD DISSIPdata 1"                                    
            print >>fs,'DISSIP 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs, DissipationEuler[i, j, k, 0], DissipationEuler[i, j, k, 1], DissipationEuler[i, j, k, 2], DissipationEuler[i, j, k, 3]                                                                                             
            print >>fs, "FIELD MASSdata 1"                                    
            print >>fs,'Mass 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,  massEuler[i, j, k, 0], massEuler[i, j, k, 1], massEuler[i, j, k, 2], massEuler[i, j, k, 3]                                                                                                                                            
            print >>fs, "FIELD EPGdata 1"                                    
            print >>fs,'EPotG 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,  epgEuler[i, j, k, 0], epgEuler[i, j, k, 1], epgEuler[i, j, k, 2], epgEuler[i, j, k, 3]                                               
            print >>fs, "FIELD EKindata 1"                                    
            print >>fs,'EKin 4 %d double'  % totalNumber
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print >>fs,  ekEuler[i, j, k, 0], ekEuler[i, j, k, 1], ekEuler[i, j, k, 2], ekEuler[i, j, k, 3]                                                                       
            fs.close()
            snap = self.read_snapshot(f)
            nfile+=1      
## --------------------------------------------------------------------
## one snapshot

class Snap:
  pass
