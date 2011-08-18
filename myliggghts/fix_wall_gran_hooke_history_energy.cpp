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
#include "fix_wall_gran_hooke_history_energy.h"
#include "pair_gran_hooke_history_energy.h"
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
#include "update.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define SMALL 1e-12

/* ---------------------------------------------------------------------- */

FixWallGranHookeHistoryEnergy::FixWallGranHookeHistoryEnergy(LAMMPS *lmp, int narg, char **arg) :
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

FixWallGranHookeHistoryEnergy::~FixWallGranHookeHistoryEnergy()
{
}

/* ---------------------------------------------------------------------- */
void FixWallGranHookeHistoryEnergy::init_substyle()
{
  //get material properties
  Yeff = ((PairGranHookeHistoryEnergy*)pairgran)->Yeff;
  Geff = ((PairGranHookeHistoryEnergy*)pairgran)->Geff;
  Kappa = ((PairGranHookeHistoryEnergy*)pairgran)->Kappa;
  betaeff = ((PairGranHookeHistoryEnergy*)pairgran)->betaeff;
  veff = ((PairGranHookeHistoryEnergy*)pairgran)->veff;
  cohEnergyDens = ((PairGranHookeHistoryEnergy*)pairgran)->cohEnergyDens;
  coeffRestLog = ((PairGranHookeHistoryEnergy*)pairgran)->coeffRestLog;
  coeffFrict = ((PairGranHookeHistoryEnergy*)pairgran)->coeffFrict;
  coeffRollFrict = ((PairGranHookeHistoryEnergy*)pairgran)->coeffRollFrict;
  charVel = ((PairGranHookeHistoryEnergy*)pairgran)->charVel;
  CPEn = ((PairGranHookeHistoryEnergy*)pairgran)->CPEn;
  CDEn = ((PairGranHookeHistoryEnergy*)pairgran)->CDEn;
  CPEt = ((PairGranHookeHistoryEnergy*)pairgran)->CPEt;
  CDEVt = ((PairGranHookeHistoryEnergy*)pairgran)->CDEVt;
  CDEFt = ((PairGranHookeHistoryEnergy*)pairgran)->CDEFt;
  CTFW = ((PairGranHookeHistoryEnergy*)pairgran)->CTFW;
  DEH = ((PairGranHookeHistoryEnergy*)pairgran)->DEH;


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
void FixWallGranHookeHistoryEnergy::reset_contact(int i,double *c_history)
{
	double &CDEnij =  c_history[3]; // this is the collision dissipated energy normal component between i and wall
	double &CDEVtij = c_history[4]; // this is the collision dissipated energy tangential component between i and wall
	double &CDEFtij = c_history[5]; // this is the collision dissipated energy tangential component between i and wall
	double &CTFWij =  c_history[6]; // this is the tangential force work term between i and wall
	DEH[i]+=(CDEnij+CDEVtij+CDEFtij+CTFWij); // The historic dissipated energy for this particle has to sum the corresponding energies for this contact that have just finished.
}


void FixWallGranHookeHistoryEnergy::updatePtrs()
{
	  CPEn = ((PairGranHookeHistoryEnergy*)pairgran)->CPEn;
	  CDEn = ((PairGranHookeHistoryEnergy*)pairgran)->CDEn;
	  CPEt = ((PairGranHookeHistoryEnergy*)pairgran)->CPEt;
	  CDEVt = ((PairGranHookeHistoryEnergy*)pairgran)->CDEVt;
	  CDEFt = ((PairGranHookeHistoryEnergy*)pairgran)->CDEFt;
	  CTFW = ((PairGranHookeHistoryEnergy*)pairgran)->CTFW;
	  DEH = ((PairGranHookeHistoryEnergy*)pairgran)->DEH;
	  if(fpgIKE) IKE=fpgIKE->values;
}

void FixWallGranHookeHistoryEnergy::compute_force(int ip, double deltan, double rsq, double meff_wall, double dx, double dy, double dz,double *vwall,double *c_history, double area_ratio)
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
  double dTx,dTy,dTz;
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

  updatePtrs();

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

  dTx=vtr1*dt;
  dTy=vtr2*dt;
  dTz=vtr3*dt;

  double &CDEnij=  c_history[3];
  double &CDEVtij= c_history[4];
  double &CDEFtij= c_history[5];
  double &CTFWij=  c_history[6];

  // well it's not clear that the pair = particle+wall can rotate
  // nevertheless, it will not harm.
  rsht = c_history[0]*dx + c_history[1]*dy + c_history[2]*dz;
  rsht = rsht*rsqinv;
  c_history[0] -= rsht*dx;
  c_history[1] -= rsht*dy;
  c_history[2] -= rsht*dz;

  double delta0x = c_history[0];
  double delta0y = c_history[1];
  double delta0z = c_history[2];
  double fe0x = -kt*c_history[0];
  double fe0y = -kt*c_history[1];
  double fe0z = -kt*c_history[2];
  double delta02 = (c_history[0]*c_history[0] + c_history[1]*c_history[1] +  c_history[2]*c_history[2]);
  double fe0    = sqrt(fe0x*fe0x+fe0y*fe0y+fe0z*fe0z);

  double dfex = - kt*dTx;
  double dfvx = - gammat*vtr1;
  double dfey = - kt*dTy;
  double dfvy = - gammat*vtr2;
  double dfez = - kt*dTz;
  double dfvz = - gammat*vtr3;

  double fe1x = fe0x+dfex;
  double fe1y = fe0y+dfey;
  double fe1z = fe0z+dfez;

  fs1 = fe0x+dfex+dfvx;
  fs2 = fe0y+dfey+dfvy;
  fs3 = fe0z+dfez+dfvz;

  fs = sqrt( fs1*fs1+fs2*fs2+fs3*fs3 );
  double fe = sqrt( fe1x*fe1x+fe1y*fe1y+fe1z*fe1z );
  fn = xmu * fabs(ccel*r);
  double fcomp = 0;
  double dfx = 0;
  double dfy = 0;
  double dfz = 0;

  if (fe>fs){
  	fcomp = fe;
  	dfx = dfex;
	dfy = dfey;
	dfz = dfez;
  }else{
  	fcomp = fs;
  	dfx = dfex+dfvx;
	dfy = dfey+dfvy;
	dfz = dfez+dfvz;
  }


//**************************************************************************

  if (fcomp > fn) {
  	if (fe0 <= fn){
  		double df2 = (dfx*dfx+dfy*dfy+dfz*dfz);
  		double lambda = (fe0x*dfx+fe0y*dfy+fe0z*dfz)/df2;
  		lambda = -lambda+sqrt(lambda*lambda+(fn*fn-fe0*fe0)/df2 );
  		if ( (lambda<0) || (lambda>1) ) error->all("Illegal value of lambda");
  		if(shearupdate){
  	  		fs1 = fe0x + lambda*dfx;
  	  		fs2 = fe0y + lambda*dfy;
  	  		fs3 = fe0z + lambda*dfz;
  	        fs = sqrt( fs1*fs1+fs2*fs2+fs3*fs3 );
  			c_history[0]+=lambda*dTx;
  			c_history[1]+=lambda*dTy;
  			c_history[2]+=lambda*dTz;
  	  		fe1x = fe0x+lambda*dfex;
  	  		fe1y = fe0y+lambda*dfey;
  	  		fe1z = fe0z+lambda*dfez;
  	  		myWorkT = -lambda*(fe1x*dTx+fe1y*dTy+fe1z*dTz);
  	  		myEdisTV= -lambda*lambda*(dfvx*dTx+dfvy*dTy+dfvz*dTz);
  	  		myEdisTF= -(1-lambda)*(fs1*dTx+fs2*dTy+fs3*dTz);
  		}
  	}else{
  		double beta = fn/fe0;
  		if(shearupdate){
  			c_history[0] *= beta;
  			c_history[1] *= beta;
  			c_history[2] *= beta;
  	  		fs1 = beta*fe0x;
  	  		fs2 = beta*fe0y;
  	  		fs3 = beta*fe0z;
  	  		myEdisTV = 0.0;
  	  		myWorkT  = -beta*(1-beta)*kt*delta02;
  	  		myEdisTF = ( (1-beta)*delta02 + (delta0x*dTx+delta0y*dTy+delta0z*dTz) )*kt*beta;
  		}
  	}
  } else {
  	if(shearupdate){
  		c_history[0] += dTx;
  		c_history[1] += dTy;
  		c_history[2] += dTz;
  		myWorkT = -(fe1x*dTx + fe1y*dTy + fe1z*dTz);
  		myEdisTV = (vtr1*vtr1+vtr2*vtr2+vtr3*vtr3)*dt*gammat;
  		myEdisTF=0.0;
  	}
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
   myEpotT = 0.0;//epK*(c_history[0]*c_history[0]+c_history[1]*c_history[1]+c_history[2]*c_history[2])*kt;

   IKE[0] += (fx*vwall[0]+fy*vwall[1]+fz*vwall[2])*dt;
//   printf("dx = %f dy = %f dz = %f\n",dx,dy,dz);
//   printf("vwx = %g vwy = %g vwz = %g IKE = %f\n",vwall[0],vwall[1],vwall[2],IKE[0]);
//   myWorkT -=(myEdisTF+myEdisTV);
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

void FixWallGranHookeHistoryEnergy::addHeatFlux(int ip, double rsq, double area_ratio)
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

inline void FixWallGranHookeHistoryEnergy::addCohesionForce(int &ip, double &r, double &Fn_coh,double area_ratio)
{
    //r is the distance between the sphere center and wall
    double Acont = (reff_wall*reff_wall-r*r)*M_PI; //contact area sphere-wall
    Fn_coh=cohEnergyDens[itype][atom_type_wall]*Acont*area_ratio;
}

/* ---------------------------------------------------------------------- */

inline void FixWallGranHookeHistoryEnergy::deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &epK)
{
    double sqrtval = sqrt(reff_wall);

    kn=16./15.*sqrtval*Yeff[itype][atom_type_wall]*pow(15.*meff_wall*charVel*charVel/(16.*sqrtval*Yeff[itype][atom_type_wall]),0.2);
    gamman=sqrt(4.*meff_wall*kn/(1.+(M_PI/coeffRestLog[itype][atom_type_wall])*(M_PI/coeffRestLog[itype][atom_type_wall])));
    switch(constflag){
    case 0:
    	kt=Kappa[itype][atom_type_wall]*kn;
    	gammat=gamman;
    	break;
    case 1:
    	kt = 2./7.*kn;
    	gammat=2./7.*gamman;
    	break;    
    case 2:
    	kt=kn;
    	gammat=gamman;
    	break;
    case 3:
    	kt=2/3*kn;
    	gammat=gamman;
    	break;
    }
    xmu=coeffFrict[itype][atom_type_wall];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall];
    epK=0.5; // this is because the integration from x
    if (dampflag == 0) gammat = 0.0;
    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE
