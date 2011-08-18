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

#include "stdlib.h"
#include "string.h"
#include "compute_grav_pe_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGPEAtom::ComputeGPEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal compute gpe/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  gpe = NULL;
  if ( (strcmp("x",arg[3])!=0 && strcmp("y",arg[3])!=0 && strcmp("z",arg[3])!=0) ) error->all("Illegal compute gpe/atom option (x/y/z)");
  char my_opcion=arg[3][0];
  switch (my_opcion)
  {
	  case 'x':
		  gpe_ncoord=0;
		  break;
	  case 'y':
		  gpe_ncoord=1;
		  break;
	  case 'z':
		  gpe_ncoord=2;
		  break;
  }
  gacc = lmp->force->numeric(arg[4]);
//  fprintf(screen,"\n\n***\n gpe_ncoord = %u\t gacc = %g\n",gpe_ncoord,gacc);
}

/* ---------------------------------------------------------------------- */

ComputeGPEAtom::~ComputeGPEAtom()
{
  memory->sfree(gpe);
}

/* ---------------------------------------------------------------------- */

void ComputeGPEAtom::init()
{
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"gpe/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute gpe/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeGPEAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(gpe);
    nmax = atom->nmax;
    gpe = (double *) memory->smalloc(nmax*sizeof(double),"gpe/atom:gpe");
    vector_atom = gpe;
  }

  // compute kinetic energy for each atom in group

  double mvv2e = force->mvv2e;
  double **x = atom->x;
//  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  if (rmass)
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
    	  gpe[i] =mvv2e * rmass[i] *gacc* x[i][gpe_ncoord];
      } else gpe[i] = 0.0;
    }

  else
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
    	  gpe[i] = mvv2e * mass[type[i]]*x[i][gpe_ncoord];
      } else gpe[i] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeGPEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
