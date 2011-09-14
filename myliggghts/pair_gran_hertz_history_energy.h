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

#ifdef PAIR_CLASS

PairStyle(gran/hertz/history/energy,PairGranHertzHistoryEnergy)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_ENERGY_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_ENERGY_H

#include "pair_gran.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryEnergy : public PairGran {

 friend class FixWallGranHertzHistoryEnergy;
 friend class FixCheckTimestepGran;

 public:

   PairGranHertzHistoryEnergy(class LAMMPS *);
  ~PairGranHertzHistoryEnergy();
  virtual void compute(int, int,int);
  virtual void settings(int, char **);
  virtual void init_substyle(); 
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

  void updatePtrs(); // Add update pointers

 protected:

  virtual void history_args(char**);
  virtual void allocate_properties(int);

  class FixPropertyGlobal* Y1; //Youngs Modulus
  class FixPropertyGlobal* v1; //Poisson's ratio
  class FixPropertyGlobal* cohEnergyDens1; //Cohesion energy density

  class FixPropertyGlobal* coeffRest1; //coefficient of restitution
  class FixPropertyGlobal* coeffFrict1; //coefficient of (static) friction
  class FixPropertyGlobal* coeffRollFrict1; //characteristic velocity needed for Linear Spring Model

  class FixPropertyGlobal* charVel1; //characteristic velocity needed for Linear Spring Model

  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,charVel,**coeffRollFrict,**Kappa;

  class FixPropertyPerAtom* fppaCPEn;
  class FixPropertyPerAtom* fppaCDEn;
  class FixPropertyPerAtom* fppaCDEVt;
  class FixPropertyPerAtom* fppaCDEFt;
  class FixPropertyPerAtom* fppaCTFW;
  class FixPropertyPerAtom* fppaDEH;
  double *CPEn, *CDEn, *CDEVt, *CDEFt, *CTFW, *DEH;

  virtual void deriveContactModelParams(int &, int &,double &, double &, double &,double &, double &, double &, double &,double &,double &);
  virtual void addCohesionForce(int &, int &,double &,double &);

  int cohesionflag; 
  int dampflag,rollingflag,constflag;
};

}

#endif
#endif
