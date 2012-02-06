import sys, commands, re, glob, types, subprocess 
from os import popen
from math import *             # any function could be used by set()
#import Numeric
from numpy import *
#from evtk.hl import gridToVTK 

#try: from DEFAULTS import PIZZA_GUNZIP
#except: PIZZA_GUNZIP = "gunzip"

# Class definition

class dumpPCOM:

    def __init__(self,*args):
        if len(args) == 0:
            raise StandardError,"dump file name not specified"
        #self.snaps = []
        #self.nsnaps =  0
        self.names = {}
        #self.atype = "type"
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
        self.x = self.names["x"]
        self.y = self.names["y"]
        self.z = self.names["z"]
        self.fx= self.names["fx"]
        self.fy= self.names["fy"]
        self.fz= self.names["fz"]
        self.vx= self.names["vx"]
        self.vy= self.names["vy"]
        self.vz= self.names["vz"]
        self.omegax=self.names["omegax"]
        self.omegay=self.names["omegay"]
        self.omegaz=self.names["omegaz"]
        self.r = self.names["radius"]    
        self.type = self.names["type"]
      
        try:
            self.cpen=self.names["f_CPEn"]
            self.cpen_flag=1
        except:
            self.cpen_flag=0 
        try:
            self.cden=self.names["f_CDEn"]
            self.cden_flag=1
        except:
            self.cden_flag=0          
        try:
            self.cdetv=self.names["f_CDEVt"]
            self.cdetv_flag=1
        except:
            self.cdetv_flag=0          
        try:
            self.cdetf=self.names["f_CDEFt"]
            self.cdetf_flag=1
        except:
            self.cdetf_flag=0                   
        try:
            self.ctfw=self.names["f_CTFW"]
            self.ctfw_flag=1
        except:
            self.ctfw_flag=0                  
        try:
            self.deh=self.names["f_DEH"]
            self.deh_flag=1
        except:
            self.deh_flag=0                       
        try:
            self.mass=self.names["mass"]
            self.mass_flag=1
        except:
            self.mass_flag=0          
        try:
            self.ePGp=self.names["v_ePGp"]
            self.ePGp_flag=1
        except:
            self.ePGp_flag=0       
        try:
            self.eKinLp=self.names["c_eKinLp"]
            self.eKinLp_flag=1
        except:
            self.eKinLp_flag=0              
        try:
            self.eKinRp=self.names["c_eKinRp"]
            self.eKinRp_flag=1
        except:
            self.eKinRp_flag=0   

        try:
            self.coordNumber=self.names["c_cn"]
            self.coordNumber_flag=1
        except:
            self.coordNumber_flag=0                     
        try:
            self.f_couple_x=self.names["f_dragforce[1]"]
            self.f_couple_y=self.names["f_dragforce[2]"]
            self.f_couple_z=self.names["f_dragforce[3]"]
            self.f_couple_flag=1
        except:
            self.f_couple_flag=0
                             
        n=0
        while snap:
            vtkFile= root+str(n).zfill(10)+".vtk"            
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
            self.writeOneVtk(fs,atoms)                    
            fs.close()            
            snap = self.read_snapshot(f)
            n+=1
    def writeAP_VTK(self, *rootF):
   
        if len(*rootF)==0:root="tmp"
        else: root=rootF[0]
      
        file = self.flist
        f=open(file)
      
        snap = self.read_snapshot(f)
        self.x = self.names["x"]
        self.y = self.names["y"]
        self.z = self.names["z"]
        self.fx= self.names["fx"]
        self.fy= self.names["fy"]
        self.fz= self.names["fz"]
        self.vx= self.names["vx"]
        self.vy= self.names["vy"]
        self.vz= self.names["vz"]
        self.omegax=self.names["omegax"]
        self.omegay=self.names["omegay"]
        self.omegaz=self.names["omegaz"]
        self.r = self.names["radius"]    
        self.type = self.names["type"]
      
        try:
            self.cpen=self.names["f_CPEn"]
            self.cpen_flag=1
        except:
            self.cpen_flag=0 
        try:
            self.cden=self.names["f_CDEn"]
            self.cden_flag=1
        except:
            self.cden_flag=0          
        try:
            self.cdetv=self.names["f_CDEVt"]
            self.cdetv_flag=1
        except:
            self.cdetv_flag=0          
        try:
            self.cdetf=self.names["f_CDEFt"]
            self.cdetf_flag=1
        except:
            self.cdetf_flag=0                   
        try:
            self.ctfw=self.names["f_CTFW"]
            self.ctfw_flag=1
        except:
            self.ctfw_flag=0                  
        try:
            self.deh=self.names["f_DEH"]
            self.deh_flag=1
        except:
            self.deh_flag=0                       
        try:
            self.mass=self.names["mass"]
            self.mass_flag=1
        except:
            self.mass_flag=0          
        try:
            self.ePGp=self.names["v_ePGp"]
            self.ePGp_flag=1
        except:
            self.ePGp_flag=0       
        try:
            self.eKinLp=self.names["c_eKinLp"]
            self.eKinLp_flag=1
        except:
            self.eKinLp_flag=0              
        try:
            self.eKinRp=self.names["c_eKinRp"]
            self.eKinRp_flag=1
        except:
            self.eKinRp_flag=0   

        try:
            self.coordNumber=self.names["c_cn"]
            self.coordNumber_flag=1
        except:
            self.coordNumber_flag=0                     
        try:
            self.f_couple_x=self.names["f_dragforce[1]"]
            self.f_couple_y=self.names["f_dragforce[2]"]
            self.f_couple_z=self.names["f_dragforce[3]"]
            self.f_couple_flag=1
        except:
            self.f_couple_flag=0
                             
        n=0
        while snap:
            vtkFileActive = root+"_active"+str(n).zfill(10)+".vtk"
            vtkFilePasive = root+"_pasive"+str(n).zfill(10)+".vtk"
            vtkFile_bb = root+str(n).zfill(10)+"_boundingBox.vtk"
            print vtkFileActive, vtkFilePasive, vtkFile_bb
            time = snap.time
            xlo=snap.xlo
            xhi=snap.xhi
            ylo=snap.ylo
            yhi=snap.yhi
            zlo=snap.zlo
            zhi=snap.zhi
            self.boundingBox(vtkFile_bb,xlo,xhi,ylo,yhi,zlo,zhi)
            atoms=snap.atoms
            atomsA = []
            atomsP = []
            for atom in atoms:
                cNumber = int(atom[self.coordNumber])
                test_activo = False
                rad_xz = math.sqrt(atom[self.x]*atom[self.x]+atom[self.z]*atom[self.z])
                try:
                    alfa = math.atan(atom[self.vz]/atom[self.vx])
                except:
                    alfa=0
                                                        
                test_activo = not((alfa<1.0) and (atom[self.z]<0.0)) and (cNumber==0) and (rad_xz<0.87)
                if test_activo:
                    atomsA.append(atom)
                else:
                    atomsP.append(atom)                        
            print 'Natoms = ', len(atoms)
            fs = open(vtkFileActive,"w")  
            self.writeOneVtk(fs,atomsA)                    
            fs.close()
            fs = open(vtkFilePasive,"w")  
            self.writeOneVtk(fs,atomsP)                    
            fs.close()          
            snap = self.read_snapshot(f)
            n+=1            
    def writeOneVtk(self, fs, atoms):

            print >>fs,"# vtk DataFile Version 2.0"
            print >>fs,"Generated by pizza.py"
            print >>fs,"ASCII"
            print >>fs,"DATASET POLYDATA"
            print >>fs,"POINTS %d float" % len(atoms)
            for atom in atoms:
                print >>fs,atom[self.x],atom[self.y],atom[self.z]  #write x,y,z  [atom[0]=id, atom[1]=type]
            print >>fs,"VERTICES",len(atoms),2*len(atoms)
            for i in xrange(len(atoms)):
                print >>fs,1,i
            print >>fs,"POINT_DATA",len(atoms)
            print >>fs,"VECTORS ","f"," float"
            for atom in atoms:
                print >>fs,atom[self.fx],atom[self.fy],atom[self.fz]
            print >>fs,"VECTORS ","v"," float"
            for atom in atoms:
                print >>fs,atom[self.vx],atom[self.vy],atom[self.vz]                        
            print >>fs,"VECTORS ","omega"," float"
            for atom in atoms:
                print >>fs,atom[self.omegax],atom[self.omegay],atom[self.omegaz]              
            print >>fs,"SCALARS atom_type int 1"
            print >>fs,"LOOKUP_TABLE default"
            for atom in atoms:  #loop all atoms
                itype = int(atom[self.type])
                print >>fs,itype,
            print >>fs,"SCALARS radius float 1"
            print >>fs,"LOOKUP_TABLE default"
            for atom in atoms:  #loop all atoms
                print >>fs,atom[self.r]  
            if self.cpen_flag:
                print >>fs,"SCALARS CPEn float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.cpen]  
            if self.cden_flag:
                print >>fs,"SCALARS CDEn float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.cden]  
            if self.cdetv_flag:
                print >>fs,"SCALARS CDEVt float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.cdetv]  
            if self.cdetf_flag:
                print >>fs,"SCALARS CDEFt float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.cdetf]                  
            if self.ctfw_flag:
                print >>fs,"SCALARS CTFW float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.ctfw]  
            if self.deh_flag:
                print >>fs,"SCALARS DEH float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.deh]                              
            if self.mass_flag:
                print >>fs,"SCALARS Mass float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.mass]                                              
            if self.ePGp_flag:
                print >>fs,"SCALARS ePGp float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.ePGp]        
            if self.eKinLp_flag:
                print >>fs,"SCALARS eKinLp float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.eKinLp]                      
            if self.eKinRp_flag:
                print >>fs,"SCALARS eKinRp float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.eKinRp]                                                      
            if self.coordNumber_flag:
                print >>fs,"SCALARS coordNumber float 1"
                print >>fs,"LOOKUP_TABLE default"
                for atom in atoms:  #loop all atoms
                    print >>fs,atom[self.coordNumber]      
                    
            if self.f_couple_flag:
                print >>fs,"VECTORS ","f_cfd"," float"
                for atom in atoms:
                    print >>fs,atom[self.f_couple_x],atom[self.f_couple_y],atom[self.f_couple_z]                    
        
    def collisionBIN(self, mycof):
       
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
      
            vn1 = delx*vnnr/rmod
            vn2 = dely*vnnr/rmod
            vn3 = delz*vnnr/rmod
        
            vt1 = vr1 - vn1;
            vt2 = vr2 - vn2;
            vt3 = vr3 - vn3;
            cr1 = r1-0.5*r1/(r1+r2)*deltan;#-0.5*deltan;
            cr2 = r2-0.5*r2/(r1+r2)*deltan;#-0.5*deltan;            
            wr1 = (cr1*atom1[omegax] + cr2*atom2[omegax]);
            wr2 = (cr1*atom1[omegay] + cr2*atom2[omegay]);
            wr3 = (cr1*atom1[omegaz] + cr2*atom2[omegaz]); 
            wrot = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
            vtr1 = vt1 - (delz*wr2-dely*wr3) /rmod;
            vtr2 = vt2 - (delx*wr3-delz*wr1) /rmod;
            vtr3 = vt3 - (dely*wr1-delx*wr2) /rmod;
            vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
            vrel = sqrt(vrel);                
            if firstTime:
                vtnorm = sqrt(vt1*vt1+vt2*vt2+vt3*vt3)
                self.t1=vt1/vtnorm
                self.t2=vt2/vtnorm
                self.t3=vt3/vtnorm
                self.vti=vtnorm
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
        self.en = fabs(vnnr/self.vni)
        self.inicio = 2.0/(1+self.en)*self.inicio/mycof
        self.final = 2.0/(1+self.en)*self.final/mycof
        self.rota = 2.0/(1+self.en)*wrot/self.vni/mycof
        if(self.vti>0):
            self.et=(vt1*self.t1+vt2*self.t2+vt3*self.t3)/self.vti
        else: self.et=1
        print self.inicio,  self.final, self.rota,  self.en, self.et
    def collisionWALL(self, myzwall, mycof):       
        firstTime=True
        file = self.flist
        f=open(file)      
        snap = self.read_snapshot(f)
        n=snap.natoms
        if n <> 1:
            raise StandardError,"Number of atoms <> 1"
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
            delx = 0
            dely = 0
            delz = atom1[z] - myzwall 
            rsq = delx*delx+dely*dely+delz*delz
            rmod=sqrt(rsq)
            r1 = atom1[r]
            deltan=r1-rmod
      
            vr1 = atom1[vx]
            vr2 = atom1[vy]
            vr3 = atom1[vz]
            vnnr = (vr1*delx+vr2*dely+vr3*delz)/rmod
      
            vn1 = delx*vnnr/rmod
            vn2 = dely*vnnr/rmod
            vn3 = delz*vnnr/rmod
        
            vt1 = vr1 - vn1;
            vt2 = vr2 - vn2;
            vt3 = vr3 - vn3;
            cr1 = r1-0.5*deltan;

            wr1 = (cr1*atom1[omegax] );
            wr2 = (cr1*atom1[omegay] );
            wr3 = (cr1*atom1[omegaz] );
            wrot = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
            vtr1 = vt1 - (delz*wr2-dely*wr3)/rmod;
            vtr2 = vt2 - (delx*wr3-delz*wr1)/rmod;
            vtr3 = vt3 - (dely*wr1-delx*wr2)/rmod;
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
        self.en = fabs(vnnr/self.vni)
        self.inicio = 2.0/(1+self.en)*self.inicio/mycof
        self.final = 2.0/(1+self.en)*self.final/mycof
        self.rota = 2.0/(1+self.en)*wrot/self.vni/mycof
        print self.inicio,  self.final ,  self.en
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
        id= self.names["id"]
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
            myID=int(atom[type]-1)
            kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
            kR = 0.5*2.0/5.0*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
            kET = kE
            kRT = kR
            cpenT = atom[cpen]            
            cdenT = atom[cden]
            cdetTV = atom[cdevt]
            cdetTF = atom[cdeft]
            ctfwT = atom[ctfw]
            dehT = atom[deh]

            EnCons = kE+kR+atom[cpen]
            EnTot = EnCons+atom[ctfw]+atom[deh]+atom[cden]+atom[cdevt]+atom[cdeft]
            EnConsT = EnCons
            EnTotT = EnTot
            print >>fs[myID], time, kE,  kR,  atom[cpen], atom[cden], atom[cdevt], atom[cdeft], atom[ctfw], atom[deh], EnCons, EnTot
            for i in xrange(1, n):
                atom=atoms[i]
                myID=int(atom[type]-1)
                kE = 0.5*atom[mass]*(atom[vx]*atom[vx]+atom[vy]*atom[vy]+atom[vz]*atom[vz])
                kR = 0.5*2.0/5.0*atom[mass]*atom[r]*atom[r]*(atom[omegax]*atom[omegax]+atom[omegay]*atom[omegay]+atom[omegaz]*atom[omegaz])
                kET += kE
                kRT += kR
                cpenT += atom[cpen]
                cdenT += atom[cden]
                cdetTV += atom[cdevt]
                cdetTF += atom[cdeft]
                ctfwT += atom[ctfw]
                dehT += atom[deh]          

                EnCons = kE+kR+atom[cpen]
                EnTot = EnCons+atom[ctfw]+atom[deh]+atom[cden]+atom[cdevt]+atom[cdeft]
                EnConsT += EnCons
                EnTotT += EnTot
                print >>fs[myID], time, kE,  kR,  atom[cpen],  atom[cden], atom[cdevt], atom[cdeft], atom[ctfw], atom[deh], EnCons, EnTot
            print >>fs[n],  time, kET, kRT, cpenT,cdenT, cdetTV, cdetTF, ctfwT, dehT, EnConsT, (EnTotT-IKE)
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
        print >>fs, '# time', 'eKIN',  'eCOL',  'ePG',  'DEH',  'IKE', 'LHS',  'ErrorA',  'ErrorR',   'kE', 'kR',  'eKIN',  'cpeN', 'cdeN',  'cdeTV',  'cdeTF',  'ctfW'
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
            cdeN = atom[cden]
            cdeTV = atom[cdevt]
            cdeTF = atom[cdeft]
            ctfW = atom[ctfw]
            eCOL= cpeN+cdeN+cdeTV+cdeTF+ctfW
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
            print >>fs, time, eKIN,  eCOL,  ePG,  DEH,  IKE,  LHS,  ErrorA,  ErrorR, kE,  kR,  eKIN,  cpeN, cdeN,  cdeTV,  cdeTF,  ctfW
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
            atoms_r = sorted(atoms,key=lambda student:student[0])
            snap.atoms = atoms_r
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
    def activeTime(self):
        
        file = self.flist
        process = subprocess.Popen(['grep', '-c', "TIMESTEP", file], stdout=subprocess.PIPE)
        stdout, stderr = process.communicate()
        timeSteps=int(stdout)                
        
        f=open(file)      
        snap = self.read_snapshot(f)
        n=snap.natoms
        id= self.names["id"]
        type = self.names["type"]
        mass=self.names["mass"]
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        vx= self.names["vx"]
        vy= self.names["vy"]
        vz= self.names["vz"]        
        fx= self.names["fx"]
        fy= self.names["fy"]
        fz= self.names["fz"]
        omegax=self.names["omegax"]
        omegay=self.names["omegay"]
        omegaz=self.names["omegaz"]
        r = self.names["radius"]
        cn = self.names["c_cn"]

        TimeArray=[]
        TimeArray=zeros((timeSteps,n+1))
        timeIndx=0
        while snap:            
            time = snap.time                  
            atomsT= snap.atoms
            print "Time = ", time
            TimeArray[timeIndx,0]=time
            particleID=1
            for atom in atomsT:                
                myID = int(atom[id])               
                Active = 0.0                 
                cNumber = int(atom[cn])
                test_activo = False
                try:
                    alfa = math.atan(atom[vz]/atom[vx])
                except:
                    alfa=0
                                                        
                test_activo = not((alfa<1.0) and (atom[z]<0.0)) and (cNumber==0)                                
                if test_activo:                    
                    Active = 1.0
                TimeArray[timeIndx,particleID]=Active
                particleID+=1
            snap = self.read_snapshot(f)
            timeIndx+=1
        Active=zeros((n))
        nombre2="activeTime.dat"
        fs=open(nombre2,"w")
        for i in range(0,timeSteps):
            print >>fs,TimeArray[i,0],
            for j in range(0,n):
                print >>fs,TimeArray[i,j+1],
                Active[j]+=TimeArray[i,j+1]
            print >>fs
        fs.close()
#        nombre2=root+"_Activetime.dat"
#        fs=open(nombre2,"w")
#         for j in range(0,n):
#            Active[j]/=timeSteps
#            print >>fs,j,Active[j]
#        fs.close()
    def flowMass(self):        
        file = self.flist
        f=open(file)      
        snap = self.read_snapshot(f)
        n=snap.natoms
        
        nombre2="flowmass.dat"
        fs=open(nombre2,"w")        
        
        id= self.names["id"]
        type = self.names["type"]
        mass=self.names["mass"]
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        vx= self.names["vx"]
        vy= self.names["vy"]
        vz= self.names["vz"]        
        fx= self.names["fx"]
        fy= self.names["fy"]
        fz= self.names["fz"]
        omegax=self.names["omegax"]
        omegay=self.names["omegay"]
        omegaz=self.names["omegaz"]
        r = self.names["radius"]
        cn = self.names["c_cn"]
        while snap:            
            ActiveFlowMass=0.0
            PassiveFlowMass=0.0
            TotalFlowMass=0.0
            time = snap.time                  
            atomsT= snap.atoms
            print "Time = ", time                    
            for atom in atomsT:                
                atomFlow=atom[mass]*atom[vy]
                TotalFlowMass+=atomFlow
                cNumber = int(atom[cn])
                test_activo = False
                rad_xz = math.sqrt(atom[x]*atom[x]+atom[z]*atom[z])
                try:
                    alfa = math.atan(atom[vz]/atom[vx])
                except:
					alfa=0
														
                test_activo = not((alfa<1.0) and (atom[z]<0.0)) and (cNumber==0) and (rad_xz<0.87)
                if test_activo:
					ActiveFlowMass+=atomFlow
                else:
					PassiveFlowMass+=atomFlow              
            print time, ActiveFlowMass, PassiveFlowMass, TotalFlowMass
            print >>fs, time, ActiveFlowMass, PassiveFlowMass, TotalFlowMass
            snap = self.read_snapshot(f)
        fs.close()

class Snap:
  pass
