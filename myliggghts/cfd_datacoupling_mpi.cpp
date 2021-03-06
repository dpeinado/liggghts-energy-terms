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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling.h"
#include "cfd_datacoupling_mpi.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

CfdDatacouplingMPI::CfdDatacouplingMPI(LAMMPS *lmp,int jarg, int narg, char **arg,FixCfdCoupling* fc) :
  CfdDatacoupling(lmp, jarg, narg, arg,fc)
{
  liggghts_is_active = false;

  this->fc = fc;

  nmax_all = 0;
  invoked_last = -1;

  allreduce_short = NULL;
  allreduce_long = NULL;

  if(comm->me == 0 && screen) fprintf(screen,"INFO: nevery as specified in LIGGGHTS is overriden by calling external program");

}

CfdDatacouplingMPI::~CfdDatacouplingMPI()
{
    memory->destroy_2d_double_array(allreduce_long);
    memory->destroy_2d_double_array(allreduce_short);
    //delete []allreduce_short;
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::pull(char *name,char *type,void *&from)
{
    
    check_grow(type,from);

    int nlocal = atom->nlocal;
    int natoms = static_cast<int>(atom->natoms);
    int len1 = -1, len2 = -1;

    void * to = find_pull_property(name,type,len1,len2);

    if (!to && nlocal)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->one("This is fatal");
    }

    if(strcmp(type,"scalar") == 0)
    {
      double **fromdouble = (double**)from;
      MPI_Allreduce(&(fromdouble[0][0]),&(allreduce_short[0][0]),natoms,MPI_DOUBLE,MPI_SUM,world);

      int m;
      double *todouble = (double*)to;
      for (int i = 0; i < natoms; i++) {
        if ((m = atom->map(i+1)) >= 0) {
          todouble[m] = allreduce_short[i][0];
        }
      }
    }

    if(strcmp(type,"vector") == 0)
    {
        double **fromdouble = (double**)from;
        MPI_Allreduce(&(fromdouble[0][0]),&(allreduce_long[0][0]),3*natoms,MPI_DOUBLE,MPI_SUM,world);
        int m;
        double **todouble = (double**)to;
        for (int i = 0; i < natoms; i++) {
          if ((m = atom->map(i+1)) >= 0) {
            todouble[m][0] = allreduce_long[i][0];
            todouble[m][1] = allreduce_long[i][1];
            todouble[m][2] = allreduce_long[i][2];
          }
        }
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::push(char *name,char *type,void *&to)
{

    check_grow(type,to);

    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int natoms = static_cast<int>(atom->natoms);

    int len1 = -1, len2 = -1;

    void * from = find_push_property(name,type,len1,len2);

    if (!from && nlocal)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s requested by calling program.\n",name);
        lmp->error->one("This is fatal");
    }

    int id;
    if(strcmp(type,"scalar") == 0)
    {
        double *fromdouble = (double*)from;
        double **todouble = (double**)to;
        for (int i = 0; i < nlocal; i++) {
            id = tag[i];
            allreduce_short[id-1][0] = fromdouble[i];
        }
        MPI_Allreduce(&(allreduce_short[0][0]),&(todouble[0][0]),natoms,MPI_DOUBLE,MPI_SUM,world);
    }
    else if(strcmp(type,"vector") == 0)
    {
        double **fromdouble = (double**)from;
        double **todouble = (double**)to;
        for (int i = 0; i < nlocal; i++) {
            id = tag[i];
            allreduce_long[id-1][0] = fromdouble[i][0];
            allreduce_long[id-1][1] = fromdouble[i][1];
            allreduce_long[id-1][2] = fromdouble[i][2];
        }
        MPI_Allreduce(&(allreduce_long[0][0]),&(todouble[0][0]),3*natoms,MPI_DOUBLE,MPI_SUM,world);
    }
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::check_grow(char *type,void *&ptr)
{
    int natoms = static_cast<int>(atom->natoms);

    if(update->ntimestep > invoked_last)
    {
        if(!atom->tag_consecutive()) error->all("MPI datatransfer model requires particles to have consecutive tags");

        growflag = false;
        while(nmax_all <= natoms)
        {
            growflag = true;
            nmax_all += DELTA;
        }

        if(growflag)
        {
            grow_arrays_allred(nmax_all);
            
        }

        invoked_last = update->ntimestep;

        reset_arrays_allred(nmax_all);
    }

    else if(growflag)
    {
        
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(int **&data, int len2,int len1,int initvalue)
{
  if(len1 == 0) len1 = atom->natoms;
  data = (int**)(lmp->memory->grow_2d_int_array(data, len1,len2, "CfdDatacouplingMPI:data"));
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

void CfdDatacouplingMPI::allocate_external(double **&data, int len2,int len1,double initvalue)
{
  if(len1 == 0) len1 = atom->natoms;
  data = (double**)(lmp->memory->grow_2d_double_array(data, len1,len2, "CfdDatacouplingMPI:data"));
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::grow_arrays_allred(int nmax)
{
  allreduce_short = (double**)(lmp->memory->grow_2d_double_array(allreduce_short, nmax,1, "fix_cfd_coupling:allreduce_short"));
  allreduce_long = (double**)(lmp->memory->grow_2d_double_array(allreduce_long, nmax,3, "fix_cfd_coupling:allreduce_long"));

  reset_arrays_allred(nmax);
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::reset_arrays_allred(int nmax)
{
  
  for (int ii = 0; ii < nmax; ii++)
      allreduce_short[ii][0] = 0.;
  for (int ii = 0; ii < nmax; ii++)
    for (int jj = 0; jj < 3; jj++)
      allreduce_long[ii][jj]=0.;
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::grow_arrays_caller(int nmax,char *type,void *&ptr)
{
   if(strcmp(type,"scalar") == 0)
   {
      ptr = (double**)(lmp->memory->grow_2d_double_array((double**)ptr, nmax,1, "fix_cfd_coupling:ptr_scalar"));
      for (int ii=0;ii<nmax;ii++)
      ((double**)ptr)[ii][0]=0.;
   }
   else if(strcmp(type,"vector") == 0)
   {
      ptr = (double**)(lmp->memory->grow_2d_double_array((double**)ptr, nmax,3, "fix_cfd_coupling:ptr_vector"));
      for (int ii=0;ii<nmax;ii++)
        for (int jj=0;jj<3;jj++)
           ((double**)ptr)[ii][jj]=0.;
   }
   else error->all("MPI datacoupling: wrong data type - must be scalar or vector") ;
}

