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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hertz_history_energyNS2.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "fix_propertyPerAtom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_propertyGlobal.h"
#include "mech_param_gran.h"
#include "compute_pair_gran_local.h"

#include <iostream>


using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHertzHistoryEnergyNS2::PairGranHertzHistoryEnergyNS2(LAMMPS *lmp) : PairGranHookeHistoryEnergyNS2(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum = 7;
    Yeff = NULL;
    Geff = NULL;
    Kappa=NULL;
    betaeff = NULL;
    veff = NULL;
    cohEnergyDens = NULL;
    coeffRestLog = NULL;
    coeffFrict = NULL;
    coeffRollFrict = NULL;
    fppaCPEn = NULL;
    fppaCDEn =  NULL;
    fppaCDEVt =  NULL;
    fppaCDEFt =  NULL;
    fppaCTFW =  NULL;
    fppaDEH =  NULL;
    CPEn = NULL;
    CDEn = NULL;
    CDEVt = NULL;
    CDEFt = NULL;
    CTFW = NULL;
    DEH = NULL;
}

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void PairGranHertzHistoryEnergyNS2::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu,double &epK)
{
    double reff=ri*rj/(ri+rj);
    double sqrtval = sqrt(reff*deltan);
    double Sn=2.*Yeff[itype][jtype]*sqrtval;
    double St=8.*Geff[itype][jtype]*sqrtval;

    kn=4./3.*Yeff[itype][jtype]*sqrtval;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);

    switch(constflag){
    case 0:
    	kt=St;
    	gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
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
    	kt=Kappa[itype][jtype]*kn;
    	gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
    	break;
    }
//    printf("\n************************************\n");
//    printf("\n kt= %f\tKappaeff = %f\n",kt,kt/kn);
//    printf("\n************************************\n");
//    printf("\n************************************\n");
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
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


/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHertzHistoryEnergyNS2::init_substyle()
{
	  int max_type = mpg->max_type();
	  allocate_properties(max_type);

	  //Get pointer to the fixes that have the material properties

	  Y1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0)]);
	  v1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0)]);

	  coeffRest1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type)]);
	  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type)]);

	  if(rollingflag)
	    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type)]);
	  if(cohesionflag)
	    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type)]);

	  //pre-calculate parameters for possible contact material combinations
	  for(int i=1;i< max_type+1; i++)
	  {
	      for(int j=1;j<max_type+1;j++)
	      {
	          double Yi=Y1->compute_vector(i-1);
	          double Yj=Y1->compute_vector(j-1);
	          double vi=v1->compute_vector(i-1);
	          double vj=v1->compute_vector(j-1);

	          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
	          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);
	          Kappa[i][j] = 2.*((1.-pow(vi,2.))/Yi  +  (1.-pow(vj,2.))/Yj)/((2.-vi)*(1.+vi)/Yi+(2.-vj)*(1.+vj)/Yj);
	          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

	          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

	          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
	          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

	          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
	          //omitting veff here

	      }
	  }

  char **fixarg = new char*[9];
  if (fppaCPEn==NULL) {
  //register Temp as property/peratom
    fixarg[0]=(char *) "CPEn";
    fixarg[1]=(char *) "all";
    fixarg[2]=(char *) "property/peratom";
    fixarg[3]=(char *) "CPEn";
    fixarg[4]=(char *) "scalar";
    fixarg[5]=(char *) "yes";
    fixarg[6]=(char *) "yes";
    fixarg[7]=(char *) "no";
    fixarg[8]=(char *) "0.0";
    modify->add_fix(9,fixarg);
    fppaCPEn=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("CPEn","property/peratom","scalar",0,0)]);
  }
  if (fppaCDEn==NULL) {
  //register Temp as property/peratom
    fixarg[0]=(char *) "CDEn";
    fixarg[1]=(char *) "all";
    fixarg[2]=(char *) "property/peratom";
    fixarg[3]=(char *) "CDEn";
    fixarg[4]=(char *) "scalar";
    fixarg[5]=(char *) "yes";
    fixarg[6]=(char *) "yes";
    fixarg[7]=(char *) "no";
    fixarg[8]=(char *) "0.0";
    modify->add_fix(9,fixarg);
    fppaCDEn=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("CDEn","property/peratom","scalar",0,0)]);
  }
  if (fppaCDEVt==NULL) {
//register Temp as property/peratom
  fixarg[0]=(char *) "CDEVt";
  fixarg[1]=(char *) "all";
  fixarg[2]=(char *) "property/peratom";
  fixarg[3]=(char *) "CDEVt";
  fixarg[4]=(char *) "scalar";
  fixarg[5]=(char *) "yes";
  fixarg[6]=(char *) "yes";
  fixarg[7]=(char *) "no";
  fixarg[8]=(char *) "0.0";
  modify->add_fix(9,fixarg);
  fppaCDEVt=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("CDEVt","property/peratom","scalar",0,0)]);
}
  if (fppaCDEFt==NULL) {
//register Temp as property/peratom
  fixarg[0]=(char *) "CDEFt";
  fixarg[1]=(char *) "all";
  fixarg[2]=(char *) "property/peratom";
  fixarg[3]=(char *) "CDEFt";
  fixarg[4]=(char *) "scalar";
  fixarg[5]=(char *) "yes";
  fixarg[6]=(char *) "yes";
  fixarg[7]=(char *) "no";
  fixarg[8]=(char *) "0.0";
  modify->add_fix(9,fixarg);
  fppaCDEFt=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("CDEFt","property/peratom","scalar",0,0)]);
}
  if (fppaCTFW==NULL) {
//register Temp as property/peratom
  fixarg[0]=(char *) "CTFW";
  fixarg[1]=(char *) "all";
  fixarg[2]=(char *) "property/peratom";
  fixarg[3]=(char *) "CTFW";
  fixarg[4]=(char *) "scalar";
  fixarg[5]=(char *) "yes";
  fixarg[6]=(char *) "yes";
  fixarg[7]=(char *) "no";
  fixarg[8]=(char *) "0.0";
  modify->add_fix(9,fixarg);
  fppaCTFW=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("CTFW","property/peratom","scalar",0,0)]);
}
  if (fppaDEH==NULL) {
 //register Temp as property/peratom
   fixarg[0]=(char *) "DEH";
   fixarg[1]=(char *) "all";
   fixarg[2]=(char *) "property/peratom";
   fixarg[3]=(char *) "DEH";
   fixarg[4]=(char *) "scalar";
   fixarg[5]=(char *) "yes";
   fixarg[6]=(char *) "yes";
   fixarg[7]=(char *) "no";
   fixarg[8]=(char *) "0.0";
   modify->add_fix(9,fixarg);
   fppaDEH=static_cast<FixPropertyPerAtom*>(modify->fix[modify->find_fix_property("DEH","property/peratom","scalar",0,0)]);
 }
  delete []fixarg;


}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */
