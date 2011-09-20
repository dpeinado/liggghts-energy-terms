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
#include "pair_gran_hertz_incremental_energyNS1.h"
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
using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHertzIncrementalEnergyNS1::PairGranHertzIncrementalEnergyNS1(LAMMPS *lmp) : PairGranHertzIncrementalEnergy(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum = 10; // 3 for previous force;  4 for CDEn, CDEVt, CDEFt, CTWF and 3 for shear
    Yeff = NULL;
    Geff = NULL;
    Kappa = NULL;
    betaeff = NULL;
    veff = NULL;
    cohEnergyDens = NULL;
    coeffRestLog = NULL;
    coeffFrict = NULL;
    coeffRollFrict = NULL;
    fppaCPEn = NULL;
    fppaCDEn = NULL;
    fppaCDEVt = NULL;
    fppaCDEFt = NULL;
    fppaCTFW = NULL;
    fppaDEH = NULL;
    CPEn = NULL;
    CDEn = NULL;
    CDEVt = NULL;
    CDEFt = NULL;
    CTFW = NULL;
    DEH = NULL;
}

/* ---------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergyNS1::compute(int eflag, int vflag, int addflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu,rmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wrmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3, meff_i,meff_j;
  double fn,fs,fs1,fs2,fs3;
  double tfs1,tfs2,tfs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dTx, dTy, dTz;

  double epK;
  double myEpotN;
  double myWorkT;
  double myEdisN;
  double myEdisTV;
  double myEdisTF;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  if (update->ntimestep > laststep) shearupdate = 1;
  else shearupdate = 0;

  updatePtrs(); // update pointers to Energy values

  fppaCPEn->do_forward_comm();
  fppaCDEn->do_forward_comm();
  fppaCDEVt->do_forward_comm();
  fppaCDEFt->do_forward_comm();
  fppaCTFW->do_forward_comm();
  fppaDEH->do_forward_comm();
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
	  i = ilist[ii];
	  CPEn[i] = 0.0; // This is equivalent to make 0.0 in force_clear()
	  CDEn[i] = 0.0; // As the contact is between i and j, and some contacts may have finished, it's needed to sum up the contacts each step
	  CDEVt[i] = 0.0; // The dissipated energy history is not cleaned. It has to accumulate all the history from beginning to the end of simulation
	  CDEFt[i] = 0.0;
	  CTFW[i] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    myEpotN = 0.;
    myEdisN = 0.;
    myEdisTV = myEdisTF = 0.;
    myWorkT = 0.;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      double mi,mj;

      if (rsq >= radsum*radsum) {
	// unset non-touching neighbors
    	if (touch[jj]){
    		touch[jj] = 0;
    		//***************************************************************
    		// Needed for introducing CDEnij to DEH[i] and DEH[j]
    	    if (rmass) {
    	      mi=rmass[i];
    	      mj=rmass[j];
    	    } else {
    	      itype = type[i];
    	      jtype = type[j];
    	      mi=mass[itype];
    	      mj=mass[jtype];
    	    }
    	    if (fr)
    	    {
    	       if(fr->body[i]>=0) double mi=fr->masstotal[fr->body[i]];
    	       if(fr->body[j]>=0) double mj=fr->masstotal[fr->body[j]];
    	    }
    	    meff=mi*mj/(mi+mj);
    	    meff_i=meff/mi;
    	    meff_j=meff/mj;
    	    if (mask[i] & freeze_group_bit) meff = mj;
    	    if (mask[j] & freeze_group_bit) meff = mi;
    	    //***************************************************************

    		shear = &allshear[dnum*jj];
    		double &CDEnij = allshear[dnum*jj+3]; // this is the collision dissipated energy normal component between i and j particles.
    		double &CDEVtij = allshear[dnum*jj+4]; // this is the collision dissipated energy tangential component between i and j particles..
    		double &CDEFtij = allshear[dnum*jj+5]; // this is the collision dissipated energy tangential component between i and j particles..
    		double &CTFWij = allshear[dnum*jj+6]; // this is the tangential force work term between i and j particles.
        	DEH[i]+=meff_i*(CDEnij+CDEVtij+CDEFtij+CTFWij); // The historic dissipated energy for this particle has to sum the corresponding energies for this contact that have just finished.
        	DEH[j]+=meff_j*(CDEnij+CDEVtij+CDEFtij+CTFWij);
    		for(int d=0; d<dnum; d++)
    			shear[d] = 0.0;
    	}
      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;
        if(!touch[jj]) touch[jj] = 1;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

	    if (rmass) {
	      mi=rmass[i];
	      mj=rmass[j];
	    } else {
	      itype = type[i];
	      jtype = type[j];
	      mi=mass[itype];
	      mj=mass[jtype];
	    }
	    if (fr)
	    {
	       if(fr->body[i]>=0) double mi=fr->masstotal[fr->body[i]];
	       if(fr->body[j]>=0) double mj=fr->masstotal[fr->body[j]];
	    }
	    meff=mi*mj/(mi+mj);
	    meff_i=meff/mi;
	    meff_j=meff/mj;
	    if (mask[i] & freeze_group_bit) meff = mj;
	    if (mask[j] & freeze_group_bit) meff = mi;

        double deltan=radsum-r;
        cri = radi;//-0.5*deltan;
        crj = radj;//-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,epK);	 //modified C.K

        damp = gamman*vnnr*rsqinv;  
        ccel = kn*(radsum-r)*rinv - damp;
        double fn_pot = kn*(radsum-r);
        //******************************************************************************************************************
        if(ccel<0) {
        	if (touch[jj]==1){
        		touch[jj] = 2;
        		myEpotN = epK*fn_pot*fn_pot/kn;
        		shear = &allshear[dnum*jj];
        		double &CDEnij = allshear[dnum*jj+3]; // this is the collision dissipated energy normal component between i and j particles.
        		double &CDEVtij = allshear[dnum*jj+4]; // this is the collision dissipated energy tangential component between i and j particles..
        		double &CDEFtij = allshear[dnum*jj+5]; // this is the collision dissipated energy tangential component between i and j particles..
        		double &CTFWij = allshear[dnum*jj+6]; // this is the tangential force work term between i and j particles.
        		DEH[i]+=meff_i*(myEpotN+CDEnij+CDEVtij+CDEFtij+CTFWij); // The historic dissipated energy for this particle has to sum the corresponding energies for this contact that have just finished.
        		DEH[j]+=meff_j*(myEpotN+CDEnij+CDEVtij+CDEFtij+CTFWij);
        		for(int d=0; d<dnum; d++)
        			shear[d] = 0.0;
            	fn_pot = 0.0;
            	ccel   = 0.0;
            	damp   = 0.0;
            	deltan = 0.0;
            	deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,epK);	 //modified C.K
        	}
        	break;
        }
        //******************************************************************************************************************
        
        if (cohesionflag) { 
            addCohesionForce(i,j,r,Fn_coh);
            ccel-=Fn_coh*rinv;
        }

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects


        double &fe0x = allshear[dnum*jj+0];
        double &fe0y = allshear[dnum*jj+1];
        double &fe0z = allshear[dnum*jj+2];
        double &CDEnij= allshear[dnum*jj+3];
        double &CDEVtij= allshear[dnum*jj+4];
        double &CDEFtij= allshear[dnum*jj+5];
        double &CTFWij= allshear[dnum*jj+6];
        shear = &allshear[dnum*jj+7];

        dTx = vtr1*dt;
    	dTy = vtr2*dt;
    	dTz = vtr3*dt;

        // The pair is rotatin as a rigid solid
    	rsht = fe0x*delx + fe0y*dely + fe0z*delz;
        rsht *= rsqinv;
        fe0x -= rsht*delx;
        fe0y -= rsht*dely;
        fe0z -= rsht*delz;
        // Now the previous tangential elastic force is in the new tangential plane

        rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
    	rsht *= rsqinv;
    	shear[0] -= rsht*delx;
    	shear[1] -= rsht*dely;
    	shear[2] -= rsht*delz;


    	double delta0x = shear[0];
        double delta0y = shear[1];
        double delta0z = shear[2];
        double delta02 = (delta0x*delta0x+delta0y*delta0y+delta0z*delta0z);
        double fe0    = sqrt(fe0x*fe0x+fe0y*fe0y+fe0z*fe0z);

    	double dfex = - kt*dTx;
    	double dfvx = - gammat*vtr1;
    	double dfey = - kt*dTy;
    	double dfvy = - gammat*vtr2;
    	double dfez = - kt*dTz;
    	double dfvz = - gammat*vtr3;

   		double fe1x = fe0x+dfex;
    	double fe1y = fe0y+dfey;
    	double fe1z = fe0z+dfez;

    	fs1 = fe0x+dfex+dfvx;
    	fs2 = fe0y+dfey+dfvy;
    	fs3 = fe0z+dfez+dfvz;

        fs = sqrt( fs1*fs1+fs2*fs2+fs3*fs3 );
        double fe = sqrt( fe1x*fe1x+fe1y*fe1y+fe1z*fe1z );
        fn = xmu * fabs(ccel*r);
        double fcomp = 0;
		double dfx = 0;
		double dfy = 0;
		double dfz = 0;

        if (fe>fs){
        	fcomp = fe;
        	dfx = dfex;
    		dfy = dfey;
    		dfz = dfez;
        }else{
        	fcomp = fs;
        	dfx = dfex+dfvx;
    		dfy = dfey+dfvy;
    		dfz = dfez+dfvz;
        }

        if (fcomp > fn) {
        	if (fe0 <= fn){
        		double df2 = (dfx*dfx+dfy*dfy+dfz*dfz);
        		double lambda = (fe0x*dfx+fe0y*dfy+fe0z*dfz)/df2;
        		lambda = -lambda+sqrt(lambda*lambda+(fn*fn-fe0*fe0)/df2 );
        		if ( (lambda<0) || (lambda>1) ) error->all("Illegal value of lambda");
                if(shearupdate){
            		fs1 = fe0x + lambda*dfx;
            		fs2 = fe0y + lambda*dfy;
            		fs3 = fe0z + lambda*dfz;
                    fs = sqrt( fs1*fs1+fs2*fs2+fs3*fs3 );
                	shear[0]+=lambda*dTx;
                	shear[1]+=lambda*dTy;
                	shear[2]+=lambda*dTz;
            		fe1x = fe0x+lambda*dfex;
            		fe1y = fe0y+lambda*dfey;
            		fe1z = fe0z+lambda*dfez;
            		myWorkT = -lambda*(fe1x*dTx+fe1y*dTy+fe1z*dTz);
            		myEdisTV= -lambda*lambda*(dfvx*dTx+dfvy*dTy+dfvz*dTz);
            		myEdisTF= -(1-lambda)*(fs1*dTx+fs2*dTy+fs3*dTz);
                }
        	}else{
        		double beta = fn/fe0;
        		if ( (beta<0) || (beta>1) ) error->all("Illegal value of beta");
        		if(shearupdate){
            		fe0x *= beta;
            		fe0y *= beta;
            		fe0z *= beta;
            		fs1 = fe0x;
            		fs2 = fe0y;
            		fs3 = fe0z;
            		myEdisTV = 0.0;
            		myWorkT  = (1-beta)*(delta0x*fe0x+delta0y*fe0y+delta0z*fe0z);
            		myEdisTF = ( -(dTx*fe0x+dTy*fe0y+dTz*fe0z) -(1-beta)*(delta0x*fe0x+delta0y*fe0y+delta0z*fe0z));
        			shear[0] = -fe0x/kt;
        		    shear[1] = -fe0y/kt;
        		    shear[2] = -fe0z/kt;
        		}
        	}
        } else {
            if(shearupdate){
        		fe0x += dfex;
        		fe0y += dfey;
        		fe0z += dfez;
            	myWorkT = -(fe0x*dTx + fe0y*dTy + fe0z*dTz);
                myEdisTV = (vtr1*vtr1+vtr2*vtr2+vtr3*vtr3)*dt*gammat;
                myEdisTF=0.0;
            	shear[0] += dTx;
            	shear[1] += dTy;
            	shear[2] += dTz;
            }
        }


        // forces & torques
        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;
        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);

        // Energy terms
        myEpotN = epK*fn_pot*fn_pot/kn; // 2/5
    	myEdisN = damp*vnnr*dt;

    	CDEnij +=  myEdisN;
    	CDEVtij += myEdisTV;
    	CDEFtij += myEdisTF;
    	CTFWij +=  myWorkT;
		if (CTFWij < 0) printf("WorkT NEGATIVO: tiempo = %u\n", update->ntimestep);
    	CPEn[i] += meff_i*myEpotN;
    	CDEn[i]+=  meff_i*CDEnij;
    	CDEVt[i]+= meff_i*CDEVtij;
    	CDEFt[i]+= meff_i*CDEFtij;
    	CTFW[i]+=  meff_i*CTFWij;

        if(rollingflag)
        {
            wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
            if (wrmag > 0.)
            {
                tor1 += rmu*kn*deltan*wr1/wrmag;
                tor2 += rmu*kn*deltan*wr2/wrmag;
                tor3 += rmu*kn*deltan*wr3/wrmag;
            }
        }
        if(addflag)
        {
        	f[i][0] += fx;
        	f[i][1] += fy;
        	f[i][2] += fz;
        	torque[i][0] -= cri*tor1;
        	torque[i][1] -= cri*tor2;
        	torque[i][2] -= cri*tor3;
        }

        if (j < nlocal) {
        	CPEn[j] += meff_j*myEpotN;
        	CDEn[j] += meff_j*CDEnij;
        	CDEVt[j]+= meff_j*CDEVtij;
        	CDEFt[j]+= meff_j*CDEFtij;
        	CTFW[j]+=  meff_j*CTFWij;
          if(addflag){
        	  f[j][0] -= fx;
        	  f[j][1] -= fy;
        	  f[j][2] -= fz;
        	  torque[j][0] -= crj*tor1;
        	  torque[j][1] -= crj*tor2;
        	  torque[j][2] -= crj*tor3;
          }
        }
        if(cpl && !addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);
        if (evflag) ev_tally_xyz(i,j,nlocal,0,
                     0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }
}
