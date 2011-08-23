/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran_hertz_incremental_energy.h"
#include "pair_gran_hertz_incremental_energy.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "fix_propertyGlobal.h"
#include "fix_propertyPerAtom.h"
#include "mech_param_gran.h"
#include "fix_rigid.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define SMALL 1e-12

/* ---------------------------------------------------------------------- */

FixWallGranHertzIncrementalEnergy::FixWallGranHertzIncrementalEnergy(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg)
{
    //args 3 and 4 are reserved for this class

    // wall/particle coefficients
    dampflag = atoi(arg[3]) & 1;
    rollingflag = atoi(arg[3]) & 2;
    cohesionflag = atoi(arg[4]) &1;
    constflag = (atoi(arg[4])&2)/2+2*(atoi(arg[4])&4)/4;
    if (cohesionflag < 0 || cohesionflag > 1 || dampflag < 0 || dampflag > 3 || constflag>3)
      error->all("Illegal fix wall/gran command");
}

/* ---------------------------------------------------------------------- */

FixWallGranHertzIncrementalEnergy::~FixWallGranHertzIncrementalEnergy()
{

}

/* ---------------------------------------------------------------------- */

void FixWallGranHertzIncrementalEnergy::init_substyle()
{
  //get material properties
  Yeff = ((PairGranHertzIncrementalEnergy*)pairgran)->Yeff;
  Geff = ((PairGranHertzIncrementalEnergy*)pairgran)->Geff;
  Kappa = ((PairGranHertzIncrementalEnergy*)pairgran)->Kappa;
  betaeff = ((PairGranHertzIncrementalEnergy*)pairgran)->betaeff;
  veff = ((PairGranHertzIncrementalEnergy*)pairgran)->veff;
  cohEnergyDens = ((PairGranHertzIncrementalEnergy*)pairgran)->cohEnergyDens;
  coeffRestLog = ((PairGranHertzIncrementalEnergy*)pairgran)->coeffRestLog;
  coeffFrict = ((PairGranHertzIncrementalEnergy*)pairgran)->coeffFrict;
  coeffRollFrict = ((PairGranHertzIncrementalEnergy*)pairgran)->coeffRollFrict;

  CPEn = ((PairGranHertzIncrementalEnergy*)pairgran)->CPEn;
  CDEn = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEn;
  CPEt = ((PairGranHertzIncrementalEnergy*)pairgran)->CPEt;
  CDEVt = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEVt;
  CDEFt = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEFt;
  CTFW = ((PairGranHertzIncrementalEnergy*)pairgran)->CTFW;
  DEH = ((PairGranHertzIncrementalEnergy*)pairgran)->DEH;


  //need to check properties for rolling friction and cohesion energy density here
  //since these models may not be active in the pair style
  
  int max_type = pairgran->mpg->max_type();
  FixPropertyGlobal *coeffRollFrict1, *cohEnergyDens1;
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type)]);
  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type)]);

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);
          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }

  if(cohesionflag) error->warning("Cohesion model should only be used with hertzian contact laws.");
}

/* ---------------------------------------------------------------------- */
void FixWallGranHertzIncrementalEnergy::reset_contact(int i,double *c_history)
{
	double &CDEnij =  c_history[3]; // this is the collision dissipated energy normal component between i and wall
	double &CDEVtij = c_history[4]; // this is the collision dissipated energy tangential component between i and wall
	double &CDEFtij = c_history[5]; // this is the collision dissipated energy tangential component between i and wall
	double &CTFWij =  c_history[6]; // this is the tangential force work term between i and wall
	DEH[i]+=(CDEnij+CDEVtij+CDEFtij+CTFWij); // The historic dissipated energy for this particle has to sum the corresponding energies for this contact that have just finished.
}
void FixWallGranHertzIncrementalEnergy::compute_force(int ip, double deltan, double rsq,double meff_wall, double dx, double dy, double dz,double *vwall,double *c_history, double area_ratio)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wrmag;
  double wr1,wr2,wr3,meff,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3;
  double shrmag,rsht,rinv,rsqinv;
  double kn, kt, gamman, gammat, xmu, rmu;
  double cri, crj;

  double *f = atom->f[ip];
  double *torque = atom->torque[ip];
  double *v = atom->v[ip];
  double *omega = atom->omega[ip];
  double radius = atom->radius[ip];
  double mass = atom->rmass[ip];
  double cr = radius - 0.5*deltan;
  double dT1,dT2,dT3;
  double epK;
  double myEpotN;
  double myEpotT;
  double myWorkT;
  double myEdisN;
  double myEdisTV;
  double myEdisTF;
  //get the parameters needed to resolve the contact
  deriveContactModelParams(ip,deltan,meff_wall,kn,kt,gamman,gammat,xmu,rmu,epK);
//  printf("\n*************************\nKN = %f\t KT = %f",kn,kt);
  CPEn = ((PairGranHertzIncrementalEnergy*)pairgran)->CPEn;
  CDEn = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEn;
  CPEt = ((PairGranHertzIncrementalEnergy*)pairgran)->CPEt;
  CDEVt = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEVt;
  CDEFt = ((PairGranHertzIncrementalEnergy*)pairgran)->CDEFt;
  CTFW = ((PairGranHertzIncrementalEnergy*)pairgran)->CTFW;
  DEH = ((PairGranHertzIncrementalEnergy*)pairgran)->DEH;
  if(fpgIKE) IKE=fpgIKE->values;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity
  // in case of wall contact, r is the contact radius

  wr1 = cr*omega[0] * rinv;
  wr2 = cr*omega[1] * rinv;
  wr3 = cr*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = gamman*vnnr*rsqinv;       
  ccel = kn*(radius-r)*rinv - damp;
  double fn_pot = kn*(radius-r);
  
  if(cohesionflag)
  {
      double Fn_coh;
      addCohesionForce(ip, r, Fn_coh,area_ratio);
      ccel-=Fn_coh*rinv;
  }

  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // shear history effects

  dT1 = vtr1*dt;
  dT2 = vtr2*dt;
  dT3 = vtr3*dt;

  double &fsix =   c_history[0];
  double &fsiy =   c_history[1];
  double &fsiz =   c_history[2];
  double &CDEnij=  c_history[3];
  double &CDEVtij= c_history[4];
  double &CDEFtij= c_history[5];
  double &CTFWij=  c_history[6];

  if(shearupdate){
	  fsix += -(kt*dT1);
	  fsiy += -(kt*dT2);
	  fsiz += -(kt*dT3);
	  // rotate shear force.
	  rsht = fsix*dx + fsiy*dy + fsiz*dz;
	  rsht *= rsqinv;
	  fsix -= rsht*dx;
	  fsiy -= rsht*dy;
	  fsiz -= rsht*dz;
    // the gammat*vtri is already in the tangential plane.
	  fs1 = fsix - gammat*vtr1;
	  fs2 = fsiy - gammat*vtr2;
	  fs3 = fsiz - gammat*vtr3;
  }
  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
  fn = xmu * fabs(ccel*r);

  myEdisTV = (vtr1*vtr1+vtr2*vtr2+vtr3*vtr3)*dt*gammat;
  myEdisTF=0.0;
  if (fs > fn) {
  	    fs1 *= fn/fs;
  	    fs2 *= fn/fs;
  	    fs3 *= fn/fs;
        fsix = fs1+gammat*vtr1;
        fsiy = fs2+gammat*vtr2;
        fsiz = fs3+gammat*vtr3;
      	myEdisTF = fabs(fs1*dT1 + fs2*dT2 + fs3*dT3);
      	myEdisTV=0.0;
  }

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  f[0] += fx;
  f[1] += fy;
  f[2] += fz;

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);
  // Energy terms
   myEpotN = epK*fn_pot*fn_pot/kn; // 0.4/2 = 2/5*1/2
   myEdisN = damp*vnnr*dt;
   myEpotT = 0.0;//0.5*(c_history[0]*c_history[0]+c_history[1]*c_history[1]+c_history[2]*c_history[2])*kt;
   myWorkT = -(fs1*dT1 + fs2*dT2 + fs3*dT3);
   IKE[0] += (fx*vwall[0]+fy*vwall[1]+fz*vwall[2])*dt;
   myWorkT -=(myEdisTF+myEdisTV);
   CDEnij +=  myEdisN;
   CDEVtij += myEdisTV;
   CDEFtij += myEdisTF;
   CTFWij +=  myWorkT;
   CPEn[ip] += myEpotN;
   CPEt[ip] += myEpotT;
   CDEn[ip]+=  CDEnij;
   CDEVt[ip]+= CDEVtij;
   CDEFt[ip]+= CDEFtij;
   CTFW[ip]+=  CTFWij;//-CDEFtij;
  if(rollingflag)
  {
	    wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
	    if (wrmag > 0.)
	    {
	        tor1 += rmu*kn*(radius-r)*wr1/wrmag;
            tor2 += rmu*kn*(radius-r)*wr2/wrmag;
            tor3 += rmu*kn*(radius-r)*wr3/wrmag;
	    }
  }

  torque[0] -= cr*tor1;
  torque[1] -= cr*tor2;
  torque[2] -= cr*tor3;
}

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

void FixWallGranHertzIncrementalEnergy::addHeatFlux(int ip, double rsq, double area_ratio)
{
    //r is the distance between the sphere center and wall
    double tcop, tcowall, hc, Acont, delta_n, r;

    r = sqrt(rsq);

    if(deltan_ratio)
    {
       delta_n = ri - r;
       delta_n *= deltan_ratio[itype-1][atom_type_wall-1];
       r = ri - delta_n;
    }

    Acont = (reff_wall*reff_wall-r*r)*M_PI*area_ratio; //contact area sphere-wall
    tcop = th_cond[itype-1]; //types start at 1, array at 0
    tcowall = th_cond[atom_type_wall-1];

    if ((fabs(tcop) < SMALL) || (fabs(tcowall) < SMALL)) hc = 0.;
    else hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(Acont);

    heatflux[ip] += (Temp_wall-Temp_p[ip]) * hc;
    
}

inline void FixWallGranHertzIncrementalEnergy::addCohesionForce(int &ip, double &r, double &Fn_coh,double area_ratio)
{
    //r is the distance between the sphere center and wall
    double Acont = (reff_wall*reff_wall-r*r)*M_PI; //contact area sphere-wall
    Fn_coh=cohEnergyDens[itype][atom_type_wall]*Acont*area_ratio;
}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHertzIncrementalEnergy::deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &epK)
{
    double sqrtval = sqrt(reff_wall*deltan);

    double Sn=2.*Yeff[itype][atom_type_wall]*sqrtval;
    double St=8.*Geff[itype][atom_type_wall]*sqrtval;

    kn=4./3.*Yeff[itype][atom_type_wall]*sqrtval;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][atom_type_wall]*sqrt(Sn*meff_wall);

    if (constflag){
    	kt=St;
    	gammat=-2.*sqrtFiveOverSix*betaeff[itype][atom_type_wall]*sqrt(St*meff_wall);
    }
    else {
    	kt = 2./7.*kn;
    	gammat=2./7.*gamman;
    }

    xmu=coeffFrict[itype][atom_type_wall];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall];
    if (dampflag == 0) gammat = 0.0;
    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    epK = 0.40;
    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE
