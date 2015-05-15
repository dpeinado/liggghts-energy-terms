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
#include "pair_gran_hertz_incremental_energy.h"
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

PairGranHertzIncrementalEnergy::PairGranHertzIncrementalEnergy(LAMMPS *lmp) : PairGran(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    /* 3 for previous tangential force elastic component;  
     * 4 for CDEn, CDEVt, CDEFt, CTWF; 
     * 1 previous normal force times xmu;
     * 1 previous deltan  
     * shear is not needed
     */
    dnum = 9;     
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

PairGranHertzIncrementalEnergy::~PairGranHertzIncrementalEnergy()
{
	//unregister energy terms as property/peratom
	  if (fppaCPEn!=NULL) modify->delete_fix("CPEn");
	  if (fppaCDEn!=NULL) modify->delete_fix("CDEn");
	  if (fppaCDEVt!=NULL) modify->delete_fix("CDEVt");
	  if (fppaCDEFt!=NULL) modify->delete_fix("CDEFt");
	  if (fppaCTFW!=NULL) modify->delete_fix("CTFW");
	  if (fppaDEH!=NULL) modify->delete_fix("DEH");
}

/* ---------------------------------------------------------------------- */
void PairGranHertzIncrementalEnergy::updatePtrs()
{

	CPEn=fppaCPEn->vector_atom;
	CDEn=fppaCDEn->vector_atom;
	CDEVt=fppaCDEVt->vector_atom;
	CDEFt=fppaCDEFt->vector_atom;
	CTFW=fppaCTFW->vector_atom;
	DEH = fppaDEH->vector_atom;

}
void PairGranHertzIncrementalEnergy::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value
    args[0]  = (char *) "fsx";
    args[1]  = (char *) "1";
    args[2]  = (char *) "fsy";
    args[3]  = (char *) "1";
    args[4]  = (char *) "fsz";
    args[5]  = (char *) "1";
    args[6] = (char *) "CDEnij";
    args[7] = (char *) "0";
    args[8] = (char *) "CDEVtij";
    args[9] = (char *) "0";
    args[10] = (char *) "CDEFtij";
    args[11] = (char *) "0";
    args[12] = (char *) "CTFWij";
    args[13] = (char *) "0";
    args[14] = (char *) "fn0";
    args[15] = (char *) "0";
}

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void PairGranHertzIncrementalEnergy::addCohesionForce(int &ip, int &jp,double &r, double &Fn_coh)
{
    //r is the distance between the sphere's centeres
    double Acont = - M_PI/4 * ( (r-ri-rj)*(r+ri-rj)*(r-ri+rj)*(r+ri+rj) )/(r*r); //contact area of the two spheres
    
    Fn_coh=cohEnergyDens[itype][jtype]*Acont;
}

/* ---------------------------------------------------------------------- */

inline void PairGranHertzIncrementalEnergy::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu, double &epK)
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
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;
    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    epK = 0.40;
    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

/* ---------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergy::compute(int eflag, int vflag, int addflag)
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
  double rsht, cri, crj;
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
    	       if(fr->body[i]>=0)  mi=fr->masstotal[fr->body[i]];
    	       if(fr->body[j]>=0)  mj=fr->masstotal[fr->body[j]];
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
	       if(fr->body[i]>=0) mi=fr->masstotal[fr->body[i]];
	       if(fr->body[j]>=0) mj=fr->masstotal[fr->body[j]];
	    }
	    meff=mi*mj/(mi+mj);
	    meff_i=meff/mi;
	    meff_j=meff/mj;
	    if (mask[i] & freeze_group_bit) meff = mj;
	    if (mask[j] & freeze_group_bit) meff = mi;

        double deltan=radsum-r;
        cri = radi-0.5*radi/(radi+radj)*deltan;
        crj = radj-0.5*radj/(radi+radj)*deltan;
        //cri = radi;//-0.5*deltan;
        //crj = radj;//-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,epK);	 //modified C.K

        damp = gamman*vnnr*rsqinv;
        ccel = kn*deltan*rinv - damp;
        double fn_pot = kn*deltan;
        //******************************************************************************************************************
        if(ccel<0) {
            /* if Fn < 0 the collision has ended. The first step after that
             * the collission disipation is transfered to the historic term
             * and everyting else is set to 0.0
             */
        	if (touch[jj]==1){
        		touch[jj] = 2;
        		myEpotN = epK*fn_pot*fn_pot/kn;
        		shear = &allshear[dnum*jj];
        		double &CDEnij = allshear[dnum*jj+3]; // this is the collision dissipated energy normal component between i and j particles.
        		double &CDEVtij = allshear[dnum*jj+4]; // this is the collision dissipated energy viscous tangential component between i and j particles..
        		double &CDEFtij = allshear[dnum*jj+5]; // this is the collision dissipated energy frictional tangential component between i and j particles..
        		double &CTFWij = allshear[dnum*jj+6]; // this is the tangential force work term between i and j particles.
        		DEH[i]+=meff_i*(myEpotN+CDEnij+CDEVtij+CDEFtij+CTFWij); // The historic dissipated energy for this particle has to sum the corresponding energies for this contact that have just finished.
        		DEH[j]+=meff_j*(myEpotN+CDEnij+CDEVtij+CDEFtij+CTFWij);
        		for(int d=0; d<dnum; d++)
        			shear[d] = 0.0;
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
        double &fn0 = allshear[dnum*jj+7];
        double &deltan0 = allshear[dnum*jj+8];

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
        // TODO falta por guardar el kt del paso anterior para usarlo en este si
        // fn<fn0 entonces uso fe0?*kt/kt(t-1) en vez de fe0?
        fn = xmu * fabs(ccel*r); // this is the total normal force (elastic + dissipation)
        double fe0    = sqrt(fe0x*fe0x+fe0y*fe0y+fe0z*fe0z);
        if(fn<fn0){
            double phi = sqrt(deltan/deltan0);
            fe0x *= phi;
            fe0y *= phi;
            fe0z *= phi;
            fe0  *= phi;
        }
    	double dfex = - kt*dTx;
    	double dfvx = - gammat*vtr1;
    	double dfey = - kt*dTy;
    	double dfvy = - gammat*vtr2;
    	double dfez = - kt*dTz;
    	double dfvz = - gammat*vtr3;

   		double fe1x = fe0x+dfex;
    	double fe1y = fe0y+dfey;
    	double fe1z = fe0z+dfez;

    	fs1 = fe1x+dfvx;
    	fs2 = fe1y+dfvy;
    	fs3 = fe1z+dfvz;

        fs = sqrt( fs1*fs1+fs2*fs2+fs3*fs3 );
        double fe = sqrt( fe1x*fe1x+fe1y*fe1y+fe1z*fe1z );
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

        if (fcomp >= fn) {
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
            		fe0x = fe0x+lambda*dfex;
            		fe0y = fe0y+lambda*dfey;
            		fe0z = fe0z+lambda*dfez;
            		myWorkT = -lambda*(fe0x*dTx+fe0y*dTy+fe0z*dTz);
            		myEdisTV= -lambda*lambda*(dfvx*dTx+dfvy*dTy+dfvz*dTz);
            		myEdisTF= -(1-lambda)*(fs1*dTx+fs2*dTy+fs3*dTz);
                }
        	}else{
        		double beta = (fn)/fe0;
        		if ( (beta<0) || (beta>1) ) error->all("Illegal value of beta");
        		if(shearupdate){
        			double dTex = ((1-beta)*fe0x)/kt;
        			double dTey = ((1-beta)*fe0y)/kt;
        			double dTez = ((1-beta)*fe0z)/kt;
        			double dTsx = dTx-dTex;
        			double dTsy = dTy-dTey;
        			double dTsz = dTz-dTez;
              		fe0x *= beta;
                	fe0y *= beta;
                	fe0z *= beta;
                	fs1 = fe0x;
                	fs2 = fe0y;
                	fs3 = fe0z;
        		    myEdisTV = 0.0;
        		    myWorkT =  -(fe0x*dTex+fe0y*dTey+fe0z*dTez);
        		    myEdisTF = -(fe0x*dTsx+fe0y*dTsy+fe0z*dTsz);
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
            }
        }
        fn0 = fn;  //Update previous value of xmu*Fn
        deltan0 = deltan; //Update previous value of deltan
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

    	CDEnij  +=  myEdisN;
    	CDEVtij += myEdisTV;
    	CDEFtij += myEdisTF;
    	CTFWij  +=  myWorkT;

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

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergy::settings(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal pair_style command");

  dampflag = force->inumeric(arg[0]) & 1;
  rollingflag = force->inumeric(arg[0]) & 2;
  cohesionflag = force->inumeric(arg[1]) & 1;;
  constflag = force->inumeric(arg[2]);
		  //(force->inumeric(arg[1])&2)/2+2*(force->inumeric(arg[1])&4)/4;
  if (dampflag < 0 || dampflag > 3 || cohesionflag < 0 || cohesionflag > 1 || constflag >3)
    error->all("Illegal pair_style command");
  if(cohesionflag && domain->dimension!=3) error->all("Cohesion model valid for 3d simulations only");
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergy::init_substyle()
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

//  charVel1=static_cast<FixPropertyGlobal*>(modify->fix[modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0)]);

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

void PairGranHertzIncrementalEnergy::allocate_properties(int size)
{
    memory->destroy_2d_double_array(Yeff);
    memory->destroy_2d_double_array(Geff);
    memory->destroy_2d_double_array(Kappa);
    memory->destroy_2d_double_array(betaeff);
    memory->destroy_2d_double_array(veff);
    memory->destroy_2d_double_array(cohEnergyDens);
    memory->destroy_2d_double_array(coeffRestLog);
    memory->destroy_2d_double_array(coeffFrict);
    memory->destroy_2d_double_array(coeffRollFrict);
    Yeff = memory->create_2d_double_array(size+1,size+1,"Yeff");
    Geff = memory->create_2d_double_array(size+1,size+1,"Geff");
    Kappa = memory->create_2d_double_array(size+1,size+1,"Kappa");
    betaeff = memory->create_2d_double_array(size+1,size+1,"betaeff");
    veff = memory->create_2d_double_array(size+1,size+1,"veff");
    cohEnergyDens = memory->create_2d_double_array(size+1,size+1,"cohEnergyDens");
    coeffRestLog = memory->create_2d_double_array(size+1,size+1,"coeffRestLog");
    coeffFrict = memory->create_2d_double_array(size+1,size+1,"coeffFrict");
    coeffRollFrict = memory->create_2d_double_array(size+1,size+1,"coeffRollFrict");
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergy::write_restart_settings(FILE *fp)
{
  fwrite(&dampflag,sizeof(int),1,fp);
  fwrite(&cohesionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts 
------------------------------------------------------------------------- */

void PairGranHertzIncrementalEnergy::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&dampflag,sizeof(int),1,fp);
    fread(&cohesionflag,sizeof(int),1,fp);
  }
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
  MPI_Bcast(&cohesionflag,1,MPI_INT,0,world);
}
