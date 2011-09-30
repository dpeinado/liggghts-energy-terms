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

PairStyle(gran/hertz/incremental/energyNS2,PairGranHertzIncrementalEnergyNS2)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_INCREMENTAL_ENERGY_NS2_H
#define LMP_PAIR_GRAN_HERTZ_INCREMENTAL_ENERGY_NS2_H

#include "pair_gran_hertz_incremental_energy.h"

namespace LAMMPS_NS {

class PairGranHertzIncrementalEnergyNS2 : public PairGranHertzIncrementalEnergy {

 friend class FixWallGranHertzIncrementalEnergy;
 friend class FixCheckTimestepGran;

 public:

  PairGranHertzIncrementalEnergyNS2(class LAMMPS *);
  virtual void compute(int, int, int);

 protected:

};

}

#endif
#endif
