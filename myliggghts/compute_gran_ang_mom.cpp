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

#include "compute_gran_ang_mom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "atom.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGranAngMom::ComputeGranAngMom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute com command");

  vector_flag = 1;
  size_vector = 3;
  extvector = 0;

  vector = new double[3];
}

/* ---------------------------------------------------------------------- */

ComputeGranAngMom::~ComputeGranAngMom()
{
}

/* ---------------------------------------------------------------------- */

void ComputeGranAngMom::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeGranAngMom::compute_vector()
{


  double **x = atom->x;
  double **v = atom->v;
  double **omega=atom->omega;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int xbox,ybox,zbox;
  double dx,dy,dz,massone,inertiaM;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double p[3];
  double lmom[3];

  invoked_vector = update->ntimestep;

  group->xcm(igroup,masstotal,Rcdg);
  group->vcm(igroup,masstotal,Vcdg);


  p[0] = p[1] = p[2] = 0.0;
  lmom[0]=lmom[1]=lmom[2]=0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i]) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      dx = (x[i][0] + xbox*xprd) - Rcdg[0];
      dy = (x[i][1] + ybox*yprd) - Rcdg[1];
      dz = (x[i][2] + zbox*zprd) - Rcdg[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      inertiaM=0.4*radius[i]*radius[i];
      p[0] += massone * (dy*(v[i][2]-Vcdg[2]) - dz*(v[i][1]-Vcdg[1])+inertiaM*omega[i][0]);
      p[1] += massone * (dz*(v[i][0]-Vcdg[0]) - dx*(v[i][2]-Vcdg[2])+inertiaM*omega[i][1]);
      p[2] += massone * (dx*(v[i][1]-Vcdg[1]) - dy*(v[i][0]-Vcdg[0])+inertiaM*omega[i][2]);
//      printf("PX: 1 = %g\t2 = %g\t3 = %g\t PX = %g\n",dy*(v[i][2]-Vcdg[2]),dz*(v[i][1]-Vcdg[1]),inertiaM*omega[i][0],p[0]);
//      printf("PY: 1 = %g\t2 = %g\t3 = %g\t PY = %g\n",dz*(v[i][0]-Vcdg[0]),dx*(v[i][2]-Vcdg[2]),inertiaM*omega[i][1],p[1]);
//      printf("PZ: 1 = %g\t2 = %g\t3 = %g\t PZ = %g\n",dx*(v[i][1]-Vcdg[1]),dy*(v[i][0]-Vcdg[0]),inertiaM*omega[i][2],p[2]);
    }
//  printf("******* PX = %g\tPY = %g\tPZ = %g\n",p[0],p[1],p[2]);
  MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
  vector=&lmom[0];
}
