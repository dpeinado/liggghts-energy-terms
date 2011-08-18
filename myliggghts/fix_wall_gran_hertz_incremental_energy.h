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

#ifdef FIX_CLASS

FixStyle(wall/gran/hertz/incremental/energy,FixWallGranHertzIncrementalEnergy)

#else

#ifndef LMP_FIX_WALL_GRAN_HERTZ_INCREMENTAL_ENERGY_H
#define LMP_FIX_WALL_GRAN_HERTZ_INCREMENTAL_ENERGY_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranHertzIncrementalEnergy : public FixWallGran {
  friend class FixTriNeighlist; 
  friend class FixConstrainMeshGran6DOF; 

 public:
  FixWallGranHertzIncrementalEnergy(class LAMMPS *, int, char **);
  ~FixWallGranHertzIncrementalEnergy();
  void updatePtrs(); // Add update pointers
 protected:

  virtual void init_substyle();
  void addHeatFlux(int, double);

  virtual void compute_force(int i,double deltan,double rsq,double meff_wall,double dx,double dy,double dz,double *vwall,double *c_history,double area_ratio);
  virtual void addHeatFlux(int i,double rsq,double area_ratio);
  virtual void addCohesionForce(int &ip, double &r, double &Fn_coh,double area_ratio);
  virtual void reset_contact(int ip,double *c_history);
  virtual void deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &epK);

  int dampflag,cohesionflag,rollingflag,constflag;
  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,**coeffRollFrict,**Kappa;
  double *CPEn, *CDEn, *CPEt, *CDEVt, *CDEFt, *CTFW, *DEH;
};

}

#endif
#endif
