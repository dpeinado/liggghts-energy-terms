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
#include "fix_wall_gran_hertz_history_energy.h"
#include "pair_gran_hertz_history_energy.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_propertyGlobal.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixWallGranHertzHistoryEnergy::FixWallGranHertzHistoryEnergy(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistoryEnergy(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */
void FixWallGranHertzHistoryEnergy::updatePtrs()
{
	  CPEn = ((PairGranHertzHistoryEnergy*)pairgran)->CPEn;
	  CDEn = ((PairGranHertzHistoryEnergy*)pairgran)->CDEn;
	  CPEt = ((PairGranHertzHistoryEnergy*)pairgran)->CPEt;
	  CDEVt = ((PairGranHertzHistoryEnergy*)pairgran)->CDEVt;
	  CDEFt = ((PairGranHertzHistoryEnergy*)pairgran)->CDEFt;
	  CTFW = ((PairGranHertzHistoryEnergy*)pairgran)->CTFW;
	  DEH = ((PairGranHertzHistoryEnergy*)pairgran)->DEH;
	  if(fpgIKE) IKE=fpgIKE->values;
}
void FixWallGranHertzHistoryEnergy::init_substyle()
{
  //get material properties
  Yeff = ((PairGranHertzHistoryEnergy*)pairgran)->Yeff;
  Geff = ((PairGranHertzHistoryEnergy*)pairgran)->Geff;
  Kappa = ((PairGranHertzHistoryEnergy*)pairgran)->Kappa;
  betaeff = ((PairGranHertzHistoryEnergy*)pairgran)->betaeff;
  veff = ((PairGranHertzHistoryEnergy*)pairgran)->veff;
  cohEnergyDens = ((PairGranHertzHistoryEnergy*)pairgran)->cohEnergyDens;
  coeffRestLog = ((PairGranHertzHistoryEnergy*)pairgran)->coeffRestLog;
  coeffFrict = ((PairGranHertzHistoryEnergy*)pairgran)->coeffFrict;
  coeffRollFrict = ((PairGranHertzHistoryEnergy*)pairgran)->coeffRollFrict;
  CPEn = ((PairGranHertzHistoryEnergy*)pairgran)->CPEn;
  CDEn = ((PairGranHertzHistoryEnergy*)pairgran)->CDEn;
  CPEt = ((PairGranHertzHistoryEnergy*)pairgran)->CPEt;
  CDEVt = ((PairGranHertzHistoryEnergy*)pairgran)->CDEVt;
  CDEFt = ((PairGranHertzHistoryEnergy*)pairgran)->CDEFt;
  CTFW = ((PairGranHertzHistoryEnergy*)pairgran)->CTFW;
  DEH = ((PairGranHertzHistoryEnergy*)pairgran)->DEH;

  //need to check properties for rolling friction and cohesion energy density here
  //since these models may not be active in the pair style
  
  FixPropertyGlobal *coeffRollFrict1, *cohEnergyDens1;
  int max_type = pairgran->mpg->max_type();
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
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void FixWallGranHertzHistoryEnergy::deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &epK)
{

    double sqrtval = sqrt(reff_wall*deltan);

    double Sn=2.*Yeff[itype][atom_type_wall]*sqrtval;
    double St=8.*Geff[itype][atom_type_wall]*sqrtval;

    kn=4./3.*Yeff[itype][atom_type_wall]*sqrtval;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][atom_type_wall]*sqrt(Sn*meff_wall);

    switch(constflag){
    case 0:
    	kt = 2./7.*kn;
    	gammat=2./7.*gamman;
    	break;
    case 1:
    	kt=St;
    	gammat=-2.*sqrtFiveOverSix*betaeff[itype][atom_type_wall]*sqrt(St*meff_wall);
    	break;
    case 2:
    	kt=kn;
    	gammat=gamman;
    	break;
    case 3:
    	kt=Kappa[itype][atom_type_wall]*kn;
    	gammat=-2.*sqrtFiveOverSix*betaeff[itype][atom_type_wall]*sqrt(St*meff_wall);
    	break;
    }

    xmu=coeffFrict[itype][atom_type_wall];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall];
    if (dampflag == 0) gammat = 0.0;
    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    epK = 0.40; // this is 2/5 because the integration from x^3/2
    return;

}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE
