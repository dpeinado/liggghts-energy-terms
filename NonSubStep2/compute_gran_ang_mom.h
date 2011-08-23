/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(gran/angmom,ComputeGranAngMom)

#else

#ifndef LMP_COMPUTE_GRAN_ANGMOM_H
#define LMP_COMPUTE_GRAN_ANGMOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGranAngMom : public Compute {
 public:
	ComputeGranAngMom(class LAMMPS *, int, char **);
  ~ComputeGranAngMom();
  void init();
  void compute_vector();

 private:
  double masstotal;
  double Rcdg[3];
  double Vcdg[3];
};

}

#endif
#endif
