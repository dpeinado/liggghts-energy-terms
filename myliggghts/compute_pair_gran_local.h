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

#ifdef COMPUTE_CLASS

ComputeStyle(pair/gran/local,ComputePairGranLocal)

#else

#ifndef LMP_COMPUTE_PAIR_GRAN_LOCAL_H
#define LMP_COMPUTE_PAIR_GRAN_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePairGranLocal : public Compute {

 public:
  ComputePairGranLocal(class LAMMPS *, int, char **);
  ~ComputePairGranLocal();
  void post_create();
  void init();
  void init_list(int, class NeighList *);
  void compute_local();
  double memory_usage();
  void add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist);
  
 private:
  int nvalues;
  int ncount;

  class PairGran *pairgran;
  int ipair;

  int posflag,idflag,fflag,tflag,hflag;

  int dnum;

  int nmax;
  double *vector;
  double **array;

  class NeighList *list;

  int count_pairs();
  void reallocate(int);
};

}

#endif
#endif
