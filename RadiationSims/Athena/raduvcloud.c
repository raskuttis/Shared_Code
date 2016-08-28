#include "copyright.h"
/*==============================================================================
 * FILE: raduvcloud.c
 *
 * Momentum injection to the ISM from star forming clouds with radiative feedback 
 * from non-ionizing UV radiation. This version initializes a spherical cloud 
 * in a box with turbulent initial conditions and follows star formation and 
 * feedback from non-ionizing UV i.e. only absorption no re-radiation in IR
 *
 * REFERENCES:
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "units.h"

#ifndef ISOTHERMAL
#error This problem generator requires --with-eos=isothermal
#endif /* ISOTHERMAL */

#define NSOLN 32768

static Real rho_small,rho_ffac,M_GMC,eps_GMC,R_GMC,kappa_IR,Psi,psi,v_final,eps_perturb;
static Real M0,tau0,tau_shell,H,Lstar,rstar,rcloud,rcrit,fcrit,tau_cloud,v_turb,surfd_out,fluxrad_out;
static Real eps_min,eps_max,rho_cloud,L_Jeans,Egrav,den_timec,rho_fdt;
static Real Lx,Ly,Lz;
static int den_counter;

/* Function prototypes for analysis and outputs in hst file */
static Real hst_dEk(const GridS *pG,const int i,const int j,const int k);
static Real hst_dEb(const GridS *pG,const int i,const int j,const int k);
static Real hst_dEb_one(const GridS *pG,const int i,const int j,const int k);
static Real hst_B3c(const GridS *pG,const int i,const int j,const int k);
static Real hst_mass_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_mass_gas(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_tot(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_gas(const GridS *pG,const int i,const int j,const int k);
static Real hst_Egrav_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_Wgrav_gas(const GridS *pG,const int i,const int j,const int k);
static Real hst_Wgrav_stars(const GridS *pG,const int i,const int j,const int k);
static Real hst_Reff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Reffone_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_rhoeffone_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Meffone_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_xeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_yeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_zeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_xxeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_yyeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_zzeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_xyeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_xzeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_yzeff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Meff_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Mass_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_MassS0_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_MassS1_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_MassS2_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_MassS3_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_MassS4_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Mr_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Efree_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_F_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_Frho_out(const GridS *pG,const int i,const int j,const int k);
static Real hst_jsrc(const GridS *pG,const int i,const int j,const int k);
static Real hst_Er(const GridS *pG,const int i,const int j,const int k);
static Real log10d(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma1(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma2(const GridS *pG,const int i,const int j,const int k);
static Real usr_Sigma3(const GridS *pG,const int i,const int j,const int k);

static Real const_absorption(const ConsS *pU);
static Real mass_dependent_luminosity(Real mass, Real age);
static Real rho(Real r);

static void diode_outflow_ix1(GridS *pG);
static void diode_outflow_ox1(GridS *pG);
static void diode_outflow_ix2(GridS *pG);
static void diode_outflow_ox2(GridS *pG);
static void diode_outflow_ix3(GridS *pG);
static void diode_outflow_ox3(GridS *pG);

/* Function prototypes for additional analysis and output e.g. PDFs and flux escaping shells */
static void intarrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, int *datarr);
static void fltarrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *datarr);
static void fltziparrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *datarr);
static void pdfarr(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *gridarr, Real *augarr, int imax, int iflag, int lflag);
static void outprocarr(const DomainS *pD, const GridS *pG, Real *locarr, Real *gridarr, Real *xcomvec, char *expr);
static Real circpdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real rmeval, Real rmaxeval, int *gblsurfpdf, Real minsd, Real delsd, int npdf, int nang, char *expr, int massflag, int radflag);
static Real anglemasspdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real rmeval, Real rmaxeval, int *gblanglepdf, int npdf);
static Real circtlimpdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real *tlims, Real rmeval, Real rmaxeval, int *gblsurfpdf, Real minsd, Real delsd, int npdf, int nang, char *expr, int massflag, int radflag);
static void fltcircpdf(const GridS *pG, Real *gblarr, Real *augarr, Real *xcomvec, Real rmeval, Real rmaxeval, Real vlim, Real *gblvpdf, int *gblsurfpdf, Real minsd, Real delsd, Real minvel, Real delvel, int npdf, int nang, int massflag);
static void fltradcircpdf(const GridS *pG, Real *augarr, Real *xcomvec, Real rmeval, Real rmaxeval, Real *gblvpdf, int *gblsurfpdf, Real minsd, Real delsd, Real minvel, Real delvel, int npdf, int nang, int massflag);
static void circfield(const GridS *pG, Real *gblarr, Real *rflux, Real *gblflux, Real *gblarea, Real *xcomvec, int nang, int nrs, int posflag);
static void loclcircfield(const GridS *pG, Real *rflux, Real *gblflux, Real *gblarea, Real *xcomvec, int nang, int nrs, int posflag, char *expr);
static Real intpgridpt(const GridS *pG, Real *gridarr, Real tx, Real ty, Real tz);
static Real intploclgridpt(const GridS *pG, Real tx, Real ty, Real tz, char *expr);
static void intupdarr(Real tval, Real minval, Real delval, int nvals, int *datarr, int updval);
static void fltupdarr(Real tval, Real minval, Real delval, int nvals, Real *datarr, Real updval);

/* Uncomment the following define to drive the flow in an impulsive manner
   as was done originally.  Restarts for this mode not yet implemented! */
/* #define IMPULSIVE_DRIVING */

/* KEEP SEMI-COLONS OUT OF THESE PRE-PROCESSOR DIRECTIVES! */
/* FFT indexing Nfast=k, Nmid=j, Nslow=i (opposite to Athena)
 * For OFST, i,j,k,nx2,nx3 reference the local grid */
#define OFST(i, j, k) ((k) + nx3*((j) + nx2*(i)))
/* KWVM: magnitude of wavenumber k in units of dkx */
#define KWVM(i, j, k) (sqrt(SQR(KCOMP(i,gis,gnx1)) +	\
                            SQR(KCOMP(j,gjs,gnx2)) +	\
                            SQR(KCOMP(k,gks,gnx3))))

/* FFTW - Variables, Plan, etc. */
/* These are made static global variables so that they need not be
   allocated AND destroyed with each call to pspect! */
static struct ath_3d_fft_plan *plan;
/* Between calls to generate(), these have unshifted, unnormalized
 * velocity perturbations. */
static ath_fft_data *fv1=NULL, *fv2=NULL, *fv3=NULL;

/* Normalized, shifted velocity perturbations */
static Real ***dv1=NULL, ***dv2=NULL, ***dv3=NULL;
/* Cutoff wavenumbers, G&O spect peak, power law spect exponent, 2 pi/L */
static Real klow,khigh,kpeak,expo,dkx;
/* Energy injection rate, last planned driving time, driving interval.
 * If not using impulsive driving, then the time quantities above are for
 * computing a new spectrum, not driving */
static Real dedt,tdrive,dtdrive;
/* Driving properties */
static int ispect,idrive;
/* Number of cells in local grid (with and without ghost zones), 
 * number of cells in global grid */
static int nx1,nx2,nx3,nx1gh,nx2gh,nx3gh,gnx1,gnx2,gnx3;
/* Starting and ending indices for local grid */
static int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
/* Starting and ending indices for global grid */
static int gis,gie,gjs,gje,gks,gke;
/* Spatial extents of local grid */
static Real x1min,x1max,x2min,x2max,x3min,x3max;
/* Length and volume elements */
static Real dx,dV;

/* Seed for random number generator */
long int rseed;
#ifdef MHD
/* beta = isothermal pressure / magnetic pressure
 * B0 = sqrt(2.0*Iso_csound2*rhobar/beta) is init magnetic field strength */
static Real beta,B0;
#endif /* MHD */
/* Initial density (will be average density throughout simulation) */
static const Real rhobar = 1.0;

/* Functions appear in this file in the same order that they appear in the
 * prototypes below */

/* Function prototypes for generating velocity perturbations */
static void pspect(ath_fft_data *ampl);
static void project();
static inline void transform();
static inline void generate();
static void perturb(GridS *pG, Real dt);

/* Function prototypes for initializing and interfacing with Athena */
static void initialize(DomainS *pD);

/* Function prototypes for Numerical Recipes functions */
static Real ran2(long int *idum);
static Real Plm(int l, int m, Real x);
static Real Ylm(int l, int m, Real theta, Real phi);
static Real factorial(int n);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int i,j,k,ii;
  int in,jn,kn;
  Real x1,x2,x3;
  Real tmp;
  Real rmin,rmax,rstart,eps,F_soln,rijk,dr,sigma,sum;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  Real dfdrcrit,A,B,C;
  Real *r_soln=NULL,*f_soln=NULL;
  Real phi,theta,Ylm_max=1.0,M_sum,M_sum_tot;

  rseed = par_geti("problem","rseed");
  if (rseed > 0.0) ath_error("[radpargrav]:  rseed must be <= 0\n");
  initialize(pD);
  tdrive = 0.0;
  
  /* Initialize gas density */
  M_sum = 0.0;
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));
        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        /* compute spherical radius */
        rijk = sqrt(SQR(x1) + SQR(x2) + SQR(x3));
        
        pG->U[k][j][i].d = MAX(rho(rijk),rho_small);
#ifdef MHD
        pG->U[k][j][i].B3c = B0;
        pG->B3i[k][j][i] = B0;
#endif
          
	/* Use passive scalars to keep track of the fluids in different regions of the cloud, since densities are same */
#if (NSCALARS > 0)
        pG->U[k][j][i].s[0] = 0.0;
        if (rijk < 0.2 * rcloud) pG->U[k][j][i].s[0] = 1.0;
        pG->U[k][j][i].s[1] = 0.0;
        if (rijk < 0.4 * rcloud && rijk >= 0.2 * rcloud) pG->U[k][j][i].s[1] = 1.0;
        pG->U[k][j][i].s[2] = 0.0;
        if (rijk < 0.6 * rcloud && rijk >= 0.4 * rcloud) pG->U[k][j][i].s[2] = 1.0;
        pG->U[k][j][i].s[3] = 0.0;
        if (rijk < 0.8 * rcloud && rijk >= 0.6 * rcloud) pG->U[k][j][i].s[3] = 1.0;
        pG->U[k][j][i].s[4] = 0.0;
        if (rijk < 1.0 * rcloud && rijk >= 0.8 * rcloud) pG->U[k][j][i].s[4] = 1.0;
#endif
          
          
      }
    }
  }
  
  /* Set the initial perturbations.  Note that we're putting in too much
   * energy this time.  This is okay since we're only interested in the
   * saturated state. */
  generate();
  perturb(pG, dtdrive);
  den_counter = 0;
  den_timec = 0.0;

  /* If decaying turbulence, no longer need the driving memory */
  if (idrive == 1) {
    ath_pout(0,"De-allocating driving memory.\n");
    
    /* Free Athena-style arrays */
    free_3d_array(dv1);
    free_3d_array(dv2);
    free_3d_array(dv3);
    
    /* Free FFTW-style arrays */
    ath_3d_fft_free(fv1);
    ath_3d_fft_free(fv2);
    ath_3d_fft_free(fv3);
  }
  
  return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *============================================================================*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pD=NULL;
  GridS *pG=NULL;
  int nl,nd;
  
  for (nl=0; nl<pM->NLevels; nl++) {
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        initialize(pD);
      }
    }
  }
  
#ifdef STAR_PARTICLE
  for (nl=0; nl<pM->NLevels; nl++) {
    for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = &(pM->Domain[nl][nd]);
        pG = pD->Grid;
        
        if (myID_Comm_world == 0) {
          starpar_printlist(0, pG);
        }
        if (pG->Gstars) {
#ifdef RADIATION
          source_exists = 1;  // PUT THIS LINE IN RESTART.C
#endif
        }
      }
    }
  }
#endif /* STAR_PARTICLE */
  
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"Sigma1")==0) return usr_Sigma1;
  if(strcmp(expr,"Sigma2")==0) return usr_Sigma2;
  if(strcmp(expr,"Sigma3")==0) return usr_Sigma3;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
    DomainS *pD=NULL;
    GridS *pG=NULL;
    int i,j,k,kk,iden,idenp,idenm,n,nc,nf,krc;
    int nl,nd,npdf=1000;
    int ntheta=100,nphi=100,nline=100, nrs=fluxrad_out, ncoms=2, nfluxes = 15, fluxind, nsdrs=surfd_out;
    Real newtime, lclrho_max, thistime, gblrho_max, di, v1, v2, v3, vmax, gblvmax, dvmax, vsq, totfloor, gblfloor, totden = 0.0, gbltotden;
    Real phimin, phimax, delphi, thetamin, thetamax, deltheta, delline, valt, svalt, cvalt, valp, svalp, cvalp, valspht, delspht;
    Real gradx, grady, gradz, tempsd, tempsdmean, tempsdouter, tx, ty, tz, tv, xd, yd, zd, nvar;
    Real temptotden, temptotdenp, temptotdenn, temptotvel, temptotvelp, temptotveln;
    Real xcomvec[3], tauout, tlims[4];
    Real xcoms, ycoms, zcoms, xcomg, ycomg, zcomg, rcomg, gblxcomg, gblycomg, gblzcomg, rcoms, rmstar, rmaxstar, totstarm, trstar, xcomt, ycomt, zcomt, rcomt;
    Real tempflux, tempden, tempxflux, tempyflux, tempzflux, temprflux, tempthetaflux;
    Real rfluxmin, rfluxmax, delrflux, minden, maxden, delden, minsd, maxsd, delsd, minvel, maxvel, delvel, delvelr, minfedd, maxfedd, delfedd;
    Real minlogfedd, maxlogfedd, dellogfedd;
    Real *allden=NULL, *gblden=NULL, *allpr=NULL, *gblpr=NULL, *rflux = NULL;
    Real *gblxflux=NULL, *gblyflux=NULL, *gblzflux=NULL, *gbldflux=NULL, *gblarea=NULL;
    Real *gblxfluxpos=NULL, *gblyfluxpos=NULL, *gblzfluxpos=NULL;
    Real *lclxflux=NULL, *lclyflux=NULL, *lclzflux=NULL, *lcldflux=NULL, *lclarea=NULL;
    Real *lclxfluxpos=NULL, *lclyfluxpos=NULL, *lclzfluxpos=NULL, *totesc=NULL;
    Real *gblvsurfpdf=NULL, *allrfluxes=NULL;
    int *gblsurfouterpdf=NULL;
    int *gblanglepdf=NULL;
    const Real rho_floor = rho_ffac*rho_small;
    char filename[24];
    char *mode, *expr;
    FILE *hp;
    int ip,jp,kp;
    int tmp,mpierr;
    
    /* If rho_floor is non-zero, enforce density/momentum floor */
    totfloor = 0.0;
    if (rho_floor) {
        for (nl=0; nl<pM->NLevels; nl++) {
            for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
                if (pM->Domain[nl][nd].Grid != NULL) {
                    pD = &(pM->Domain[nl][nd]);
                    pG = pD->Grid;
                    
                    for (k=kl; k<=ku; k++) {
                        for (j=kl; j<=ju; j++) {
                            for (i=il; i<=iu; i++) {
                                if (pG->U[k][j][i].d < rho_floor) {
                                    totfloor = totfloor + (rho_floor - pG->U[k][j][i].d);
                                    pG->U[k][j][i].d = rho_floor;
                                    pG->U[k][j][i].M1 = 0.0;
                                    pG->U[k][j][i].M2 = 0.0;
                                    pG->U[k][j][i].M3 = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Output Star Particle Data and calculate the stellar centre of mass
#ifdef STAR_PARTICLE
    StarParListS *pGstars = NULL;
    StarParS *pStar = NULL;
    
    for (nl=0; nl<pM->NLevels; nl++) {
        for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
            if (pM->Domain[nl][nd].Grid != NULL) {
                pD = &(pM->Domain[nl][nd]);
                pG = pD->Grid;
                
                // Calculate the centre of mass for star particles
                pGstars = pG->Gstars;
                totstarm = 0.0;
                xcoms = 0.0;
                ycoms = 0.0;
                zcoms = 0.0;
                rmstar = 0.0;
                rmaxstar = 0.0;
                while (pGstars) {
                    pStar = &(pGstars->starpar);
                    totstarm = totstarm + pStar->m;
                    xcoms = xcoms + pStar->x1 * pStar->m;
                    ycoms = ycoms + pStar->x2 * pStar->m;
                    zcoms = zcoms + pStar->x3 * pStar->m;
                    pGstars = pGstars->next;
                }
                if (totstarm > 0.0) {
                    xcoms = xcoms / totstarm;
                    ycoms = ycoms / totstarm;
                    zcoms = zcoms / totstarm;
                }
                rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
                
                pGstars = pG->Gstars;
                while (pGstars) {
                    pStar = &(pGstars->starpar);
                    trstar = (pStar->x1 - xcoms) * (pStar->x1 - xcoms) + (pStar->x2 - ycoms) * (pStar->x2 - ycoms) + (pStar->x3 - zcoms) * (pStar->x3 - zcoms);
                    rmstar = rmstar + trstar * pStar->m / totstarm;
                    if (trstar > rmaxstar) {
                        rmaxstar = trstar;
                    }
                    pGstars = pGstars->next;
                }
                if (totstarm > 0.0) {
                    rmstar = sqrt(rmstar);
                    rmaxstar = sqrt(rmaxstar);
                }
                
                /* Write star particle status to output file */
                if (myID_Comm_world == 0) {
                    pGstars = pG->Gstars;
                    
                    while (pGstars) {
                        pStar = &(pGstars->starpar);
                        mode = (pStar->age == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                        sprintf(filename,"../star%4.4d.dat",pStar->id);
#else
                        sprintf(filename,"star%4.4d.dat",pStar->id);
#endif
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open starpar dump file\n");
                            return;
                        }
                        fprintf(hp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
                                pG->time,pStar->m,pStar->x1,pStar->x2,pStar->x3,
                                pStar->v1,pStar->v2,pStar->v3,pStar->age,pStar->mdot,
                                pStar->merge_history);
                        fclose(hp);
                        pGstars = pGstars->next;
                    }
                }
            }
        }
    }
    
#endif /* STAR_PARTICLE */
    
    /* Only calculate and make outputs if dt has passed since the last output */
    if (pG->time >= den_timec) {
        
        // Calculate the maximum density, the gas center of mass and the total center of mass (including gas and stars)
        lclrho_max = 0.0;
        vmax = 0.0;
        dvmax = 0.0;
        xcomg = 0.0;
        ycomg = 0.0;
        zcomg = 0.0;
        totden = 0.0;
        for (nl=0; nl<pM->NLevels; nl++) {
            for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
                if (pM->Domain[nl][nd].Grid != NULL) {
                    pD = &(pM->Domain[nl][nd]);
                    pG = pD->Grid;
                    
                    for (k=ks; k<=ke; k++) {
                        // ath_pout(1,"time = %e, rho_max = %e, k = [%d, %d], processor = %d\n",pG->time, lclrho_max, k, kl, myID_Comm_world);
                        for (j=js; j<=je; j++) {
                            for (i=is; i<=ie; i++) {
                                totden = totden + pG->U[k][j][i].d;
                                cc_pos(pG,i,j,k,&tx,&ty,&tz);
                                xcomg = xcomg + pG->U[k][j][i].d * tx;
                                ycomg = ycomg + pG->U[k][j][i].d * ty;
                                zcomg = zcomg + pG->U[k][j][i].d * tz;
                                if (pG->U[k][j][i].d > lclrho_max) {
                                    lclrho_max = pG->U[k][j][i].d;
                                }
                                di = 1.0/(pG->U[k][j][i].d);
                                v1 = pG->U[k][j][i].M1*di;
                                v2 = pG->U[k][j][i].M2*di;
                                v3 = pG->U[k][j][i].M3*di;
                                vsq = SQR(v1) + SQR(v2) + SQR(v3);
                                if (vsq > vmax) {
                                    vmax = vsq;
                                    dvmax = 1.0/di;
                                }
                            }
                        }
                    }
                    
#ifdef MPI_PARALLEL
                    mpierr = MPI_Reduce(&(lclrho_max), &(gblrho_max), 1,
                                        MPI_DOUBLE, MPI_MAX, 0, pD->Comm_Domain);
                    mpierr = MPI_Reduce(&(totfloor), &(gblfloor), 1,
                                        MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                    mpierr = MPI_Reduce(&(xcomg), &(gblxcomg), 1,
                                        MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                    mpierr = MPI_Reduce(&(ycomg), &(gblycomg), 1,
                                        MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                    mpierr = MPI_Reduce(&(zcomg), &(gblzcomg), 1,
                                        MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                    mpierr = MPI_Reduce(&(totden), &(gbltotden), 1,
                                        MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif
                    gblxcomg = gblxcomg / gbltotden;
                    gblycomg = gblycomg / gbltotden;
                    gblzcomg = gblzcomg / gbltotden;
                    gbltotden = gbltotden * pG->dx1 * pG->dx2 * pG->dx3;
                    rcomg = sqrt(gblxcomg * gblxcomg + gblycomg * gblycomg + gblzcomg * gblzcomg);
                    xcomt = (gblxcomg * gbltotden + xcoms * totstarm) / (gbltotden + totstarm);
                    ycomt = (gblycomg * gbltotden + ycoms * totstarm) / (gbltotden + totstarm);
                    zcomt = (gblzcomg * gbltotden + zcoms * totstarm) / (gbltotden + totstarm);
                    rcomt = sqrt(xcomt * xcomt + ycomt * ycomt + zcomt * zcomt);
                    
                    // Output the maximum density, the density added by flooring and all centers of mass
                    if (myID_Comm_world == 0) {
                        
                        mode = (pG->time == 0.0) ? "w" : "a";
                        
#ifdef MPI_PARALLEL
                        sprintf(filename,"../denmaxall.dat");
#else
                        sprintf(filename,"denmaxall.dat");
#endif
                        
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open density dump file\n");
                            return;
                        }
                        
                        fprintf(hp,"%lf %lf\n", pG->time, gblrho_max);
                        fclose(hp);
                        
                        mode = (pG->time == 0.0) ? "w" : "a";
                        
#ifdef MPI_PARALLEL
                        sprintf(filename,"../denfloorall.dat");
#else
                        sprintf(filename,"denfloorall.dat");
#endif
                        
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open density dump file\n");
                            return;
                        }
                        
                        fprintf(hp,"%lf %lf\n", pG->time, gblfloor);
                        fclose(hp);
                        
                        // Centres of Mass
                        mode = (pG->time == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                        sprintf(filename,"../com_star.dat");
#else
                        sprintf(filename,"com_star.dat");
#endif
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open density dump file\n");
                            return;
                        }
                        
                        fprintf(hp,"\n%lf %e %e %e %e\n", pG->time, xcoms, ycoms, zcoms, rcoms);
                        fclose(hp);
                        
                        // Gas centres of mass
                        mode = (pG->time == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                        sprintf(filename,"../com_gas.dat");
#else
                        sprintf(filename,"com_gas.dat");
#endif
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open density dump file\n");
                            return;
                        }
                        
                        fprintf(hp,"\n%lf %e %e %e %e\n", pG->time, gblxcomg, gblycomg, gblzcomg, rcomg);
                        fclose(hp);
                        
                        // Total centres of mass
                        mode = (pG->time == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                        sprintf(filename,"../com_total.dat");
#else
                        sprintf(filename,"com_total.dat");
#endif
                        if((hp = fopen(filename,mode)) == NULL) {
                            ath_error("[radpargrav]: Unable to open density dump file\n");
                            return;
                        }
                        
                        fprintf(hp,"\n%lf %e %e %e %e\n", pG->time, xcomt, ycomt, zcomt, rcomt);
                        fclose(hp);
                    }
                }
            }
        }
        
        // This flag is used to indicated whether we want to output PDFs or not (1 or 0)
        if (surfd_out > 0) {
        
            // Limits for the density pdf
            minden = log10(rho_floor);
            maxden = log10(rho_cloud) + 6.0;
            delden = (maxden - minden) / ((double) (npdf) - 1.0);
            
            // Limits for the surface density pdf
            minsd = log10(rho_cloud * rcloud);
            minsd = minsd - 6.0;
            maxsd = minsd + 12.0;
            delsd = (maxsd - minsd) / ((double) (npdf) - 1.0);
            
            // Limits for the Eddington ratio pdf
            minfedd = -100.0;
            maxfedd = 100.0;
            delfedd = (maxfedd - minfedd) / ((double) (npdf) - 1.0);
            
            // Log of the same limits
            minlogfedd = -4.0;
            maxlogfedd = 4.0;
            dellogfedd = (maxlogfedd - minlogfedd) / ((double) (npdf) - 1.0);
            
            // Limits of the velocity pdf (delvelr is two-sided)
            minvel = 0.0;
            maxvel = 5.0 * v_turb;
            delvel = (maxvel - minvel) / ((double)(npdf) - 1.0);
            delvelr = (2.0 * maxvel) / ((double)(npdf) - 1.0);
            
            // Communicate the whole density array to the output processor
            if ((allden=(Real*)calloc_1d_array(gnx1*gnx2*gnx3,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for density grid\n");
            }
            if ((gblden=(Real*)calloc_1d_array(gnx1*gnx2*gnx3,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for density grid\n");
            }
            if ((allpr=(Real*)calloc_1d_array(gnx1*gnx2*gnx3,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblpr=(Real*)calloc_1d_array(gnx1*gnx2*gnx3,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblsurfouterpdf=(int*)calloc_1d_array(npdf,sizeof(int)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for surface density pdf\n");
            }
            if ((gblanglepdf=(int*)calloc_1d_array(ntheta,sizeof(int)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for surface density pdf\n");
            }
            if ((gblvsurfpdf=(Real*)calloc_1d_array(npdf,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for surface density velocity pdf\n");
            }
            
            // Calculate global arrays at each processor
            for (nl=0; nl<pM->NLevels; nl++) {
                for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
                    if (pM->Domain[nl][nd].Grid != NULL) {
                        pD = &(pM->Domain[nl][nd]);
                        pG = pD->Grid;
                        
                        // xcomvec is used to decide the origin for circumcluster PDFs
                        xcomvec[0] = xcoms;
                        xcomvec[1] = ycoms;
                        xcomvec[2] = zcoms;
                        
                        // Store full density and mod of velocity across whole grid
                        expr = "d";
                        outprocarr(pD, pG, allden, gblden, xcomvec, expr);
                        expr = "Vt";
                        outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                        
                        if (myID_Comm_world == 0) {
                            // Calculation of the surface density in the z-direction
                            for (i = 0; i <= gnx1 - 1; i++) {
                                for (j = 0; j <= gnx2 - 1; j++) {
                                    tempsd = 0.0;
                                    for (k = 0; k <= gnx3 - 1; k++) {
                                        tempsd = tempsd + gblden[i + gnx1 * j + gnx1 * gnx2 * k] * pG->dx1;
                                    }
                                    intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfouterpdf, 1);
                                }
                            }
                            
                            // Surface Densities in the z-direction
                            sprintf(filename,"sdpdfxy.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            for (i = 1; i <= npdf; i++) {
                                gblsurfouterpdf[i-1] = 0;
                            }
                            
                            // Calculation of the surface density in the x-direction
                            for (k = 0; k <= gnx3 - 1; k++) {
                                for (j = 0; j <= gnx2 - 1; j++) {
                                    tempsd = 0.0;
                                    for (i = 0; i <= gnx1 - 1; i++) {
                                        tempsd = tempsd + gblden[i + gnx1 * j + gnx1 * gnx2 * k] * pG->dx1;
                                    }
                                    intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfouterpdf, 1);
                                }
                            }
                            
                            // Surface Densities in the x-direction
                            sprintf(filename,"sdpdfyz.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            for (i = 1; i <= npdf; i++) {
                                gblsurfouterpdf[i-1] = 0;
                            }
                            
                            // Calculation of the surface density in the y-direction
                            for (i = 0; i <= gnx1 - 1; i++) {
                                for (k = 0; k <= gnx3 - 1; k++) {
                                    tempsd = 0.0;
                                    for (j = 0; j <= gnx2 - 1; j++) {
                                        tempsd = tempsd + gblden[i + gnx1 * j + gnx1 * gnx2 * k] * pG->dx1;
                                    }
                                    intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfouterpdf, 1);
                                }
                            }
                            
                            // Surface Densities in the y-direction
                            sprintf(filename,"sdpdfxz.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            for (i = 1; i <= npdf; i++) {
                                gblsurfouterpdf[i-1] = 0;
                            }
                            
                            // Circumcluster Densities Around the Stellar Centre of Mass
                            ntheta = 100;
                            tauout = circpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 0);
                            sprintf(filename,"sdpdfallcirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Surface Densities Around the Stellar Centre of Mass from the Outer Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 0);
                            sprintf(filename,"sdpdfoutercirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Surface Densities Around the Stellar Centre of Mass from the Mean Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 0);
                            sprintf(filename,"sdpdfmeancirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass
                            tauout = circpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                            sprintf(filename,"sdmasspdfallcirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass from the Outer Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                            sprintf(filename,"sdmasspdfoutercirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass from the Mean Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                            sprintf(filename,"sdmasspdfmeancirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Circumcluster densities restricted to angular bins above and below the x-y plane
                            /* for (i = 0; i < 10; i++) {
                            
                                xcomvec[0] = xcoms;
                                xcomvec[1] = ycoms;
                                xcomvec[2] = zcoms;
                                nvar = (Real)(i) / 10.0;
                                // Mass-Weighted Surface Densities Around the Stellar Centre of Mass from the Mean Stellar Radius within 45 degrees of x-y plane and y-z plane
                                tlims[0] = 3.14159 / 2.0 * nvar;
                                tlims[1] = 3.14159 / 2.0;
                                tlims[2] = 3.14159 / 2.0;
                                tlims[3] = 3.14159 / 2.0 * (2.0 - nvar);
                                // ath_pout(0,"i = %d, nvar = %e, tlims = [%e, %e, %e, %e]\n", i, nvar, tlims[0], tlims[1], tlims[2], tlims[3]);
                                tauout = circtlimpdf(pG, gblden, xcomvec, tlims, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                                sprintf(filename,"sdmpdfxyplane_%d.dat", i);
                                intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                                
                                tlims[0] = 0.0;
                                tlims[1] = 3.14159 / 2.0 * nvar;
                                tlims[2] = 3.14159 / 2.0 * (2.0 - nvar);
                                tlims[3] = 3.14159;
                                // ath_pout(0,"i = %d, nvar = %e, tlims = [%e, %e, %e, %e]\n", i, nvar, tlims[0], tlims[1], tlims[2], tlims[3]);
                                tauout = circtlimpdf(pG, gblden, xcomvec, tlims, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                                sprintf(filename,"sdmpdfyzplane_%d.dat", i);
                                intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                                
                                xcomvec[0] = xcomt;
                                xcomvec[1] = ycomt;
                                xcomvec[2] = zcomt;
                                tlims[0] = 3.14159 / 2.0 * nvar;
                                tlims[1] = 3.14159 / 2.0;
                                tlims[2] = 3.14159 / 2.0;
                                tlims[3] = 3.14159 / 2.0 * (2.0 - nvar);
                                tauout = circtlimpdf(pG, gblden, xcomvec, tlims, rmstar, 2.0 * rcloud - rcomt - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                                sprintf(filename,"sdmpdfxyplanet_%d.dat", i);
                                intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                                
                                tlims[0] = 0.0;
                                tlims[1] = 3.14159 / 2.0 * nvar;
                                tlims[2] = 3.14159 / 2.0 * (2.0 - nvar);
                                tlims[3] = 3.14159;
                                tauout = circtlimpdf(pG, gblden, xcomvec, tlims, rmstar, 2.0 * rcloud - rcomt - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 0);
                                sprintf(filename,"sdmpdfyzplanet_%d.dat", i);
                                intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                                
                            } */
                            // ath_error("Premature Exit\n");
                            
                            // PDF of mass as a function of angle above the x-y plane around the gaseous centre of mass
                            xcomvec[0] = gblxcomg;
                            xcomvec[1] = gblycomg;
                            xcomvec[2] = gblzcomg;
                            anglemasspdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcomg - pG->dx1, gblanglepdf, ntheta);
                            // For all mass
                            sprintf(filename,"angmasspdfallg.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            anglemasspdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcomg - pG->dx1, gblanglepdf, ntheta);
                            // From only the mean stellar radius outwards
                            sprintf(filename,"angmasspdfmeang.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            // Same but around the total center of mass
                            xcomvec[0] = xcomt;
                            xcomvec[1] = ycomt;
                            xcomvec[2] = zcomt;
                            anglemasspdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcomt - pG->dx1, gblanglepdf, ntheta);
                            sprintf(filename,"angmasspdfallt.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            anglemasspdf(pG, gblden, xcomvec, 10.0, 2.0 * rcloud - rcomt - pG->dx1, gblanglepdf, ntheta);
                            sprintf(filename,"angmasspdfmeant.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            // Same but around the stellar center of mass
                            xcomvec[0] = xcoms;
                            xcomvec[1] = ycoms;
                            xcomvec[2] = zcoms;
                            anglemasspdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblanglepdf, ntheta);
                            sprintf(filename,"angmasspdfall.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            anglemasspdf(pG, gblden, xcomvec, 10.0, 2.0 * rcloud - rcoms - pG->dx1, gblanglepdf, ntheta);
                            sprintf(filename,"angmasspdfmean.dat");
                            intarrdatout(pG->time, 1.0, -1.0 / (1.0 * ntheta), ntheta, filename, gblanglepdf);
                            
                            // Same, but mass-weight along the line rather than weight equally
                            /* tauout = circpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 1);
                            sprintf(filename,"mwsdpdfallcirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Surface Densities Around the Stellar Centre of Mass from the Outer Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 1);
                            sprintf(filename,"mwsdpdfoutercirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Surface Densities Around the Stellar Centre of Mass from the Mean Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 0, 1);
                            sprintf(filename,"mwsdpdfmeancirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass
                            tauout = circpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 1);
                            sprintf(filename,"mwsdmasspdfallcirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass from the Outer Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 1);
                            sprintf(filename,"mwsdmasspdfoutercirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf);
                            
                            // Mass-Weighted Surface Densities Around the Stellar Centre of Mass from the Mean Stellar Radius
                            tauout = circpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblsurfouterpdf, minsd, delsd, npdf, ntheta, "sd", 1, 1);
                            sprintf(filename,"mwsdmasspdfmeancirc.dat");
                            intarrdatout(pG->time, minsd, delsd, npdf, filename, gblsurfouterpdf); */
                            
                            // Density pdf
                            sprintf(filename,"denpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, gnx1*gnx2*gnx3, 0, 1);
                            
                            // Velocity (absolute) as a function of density
                            sprintf(filename,"vdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 1, 1);
                            sprintf(filename,"vsigmadenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 2, 1);
                            
                            // PDF of absolute value of velocity
                            // sprintf(filename,"vpdfall.dat");
                            // pdfarr(pG->time, minvel, delvel, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                            // sprintf(filename,"vmasspdfall.dat");
                            // pdfarr(pG->time, minvel, delvel, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                            
                            // PDF of absolute value of velocity over a broader range of velocities
                            sprintf(filename,"vlrpdfall.dat");
                            pdfarr(pG->time, 2.5 * minvel, 2.5 * delvel, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                            sprintf(filename,"vlrmasspdfall.dat");
                            pdfarr(pG->time, 2.5 * minvel, 2.5 * delvel, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                            
                            // PDF of absolute value of logarithm of velocity
                            // sprintf(filename,"vlogpdfall.dat");
                            // pdfarr(pG->time, -1.0, 0.004, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 1);
                            // sprintf(filename,"vlogmasspdfall.dat");
                            // pdfarr(pG->time, -1.0, 0.004, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 1);
                            
                        }
                        
                        // PDFs of outflowing velocity - Need the total array of radial momentum
                        expr = "Mr";
                        outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                        if (myID_Comm_world == 0) {
                            
                            // Radial Velocity as a function of circumcluster surface density
                            fltcircpdf(pG, gblpr, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, 0.0, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"vsdpdfallcirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the mean stellar radius
                            fltcircpdf(pG, gblpr, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, 0.0, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"vsdpdfmeancirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the outer stellar radius
                            fltcircpdf(pG, gblpr, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, 0.0, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"vsdpdfoutercirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Mean radius as a function of circumcluster surface density
                            fltradcircpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"rsdpdfallcirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the mean stellar radius
                            fltradcircpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"rsdpdfmeancirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the outer stellar radius
                            fltradcircpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 0);
                            sprintf(filename,"rsdpdfoutercirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Radial Velocity as a function of mass-weighted circumcluster surface density
                            fltcircpdf(pG, gblpr, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, -1.0e5, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"vsdmasspdfallcirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the mean stellar radius
                            fltcircpdf(pG, gblpr, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, -1.0e5, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"vsdmasspdfmeancirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the outer stellar radius
                            fltcircpdf(pG, gblpr, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, -1.0e5, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"vsdmasspdfoutercirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Mean radius as a function of mass-weighted circumcluster surface density
                            fltradcircpdf(pG, gblden, xcomvec, 0.0, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"rsdmasspdfallcirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the mean stellar radius
                            fltradcircpdf(pG, gblden, xcomvec, rmstar, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"rsdmasspdfmeancirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                            
                            // Same but for circumcluster distribution away from the outer stellar radius
                            fltradcircpdf(pG, gblden, xcomvec, rmaxstar, 2.0 * rcloud - rcoms - pG->dx1, gblvsurfpdf, gblsurfouterpdf, minsd, delsd, minvel, delvel, npdf, ntheta, 1);
                            sprintf(filename,"rsdmasspdfoutercirc.dat");
                            fltziparrdatout(pG->time, minsd, delsd, npdf, filename, gblvsurfpdf);
                        }
                        
                        // PDF of Eddington Factor
#ifdef RADIATION
                        expr = "Fedd";
                        outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                        if (myID_Comm_world == 0) {
                            sprintf(filename,"feddpdfall.dat");
                            pdfarr(pG->time, minfedd, delfedd, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                            // Mass-weighted
                            sprintf(filename,"feddmasspdfall.dat");
                            pdfarr(pG->time, minfedd, delfedd, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                            // sprintf(filename,"logfeddpdfall.dat");
                            // pdfarr(pG->time, minlogfedd, dellogfedd, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 1);
                            // sprintf(filename,"logfeddmasspdfall.dat");
                            // pdfarr(pG->time, minlogfedd, dellogfedd, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 1);
                            // Eddington Factor as a function of density
                            // PDF of Eddington factor as a function of cell density
                            sprintf(filename,"fedddenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 1, 1);
                            // PDF of Eddington factor by LOS circumcluster density
                            sprintf(filename,"feddsigmadenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 2, 1);
                        }
#endif
                        
                        // PDF of radial velocity
                        expr = "Vr";
                        outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                        if (myID_Comm_world == 0) {
                            // sprintf(filename,"vrpdfall.dat");
                            // pdfarr(pG->time, -1.0 * maxvel, delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                            // sprintf(filename,"vrmasspdfall.dat");
                            // pdfarr(pG->time, -1.0 * maxvel, delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                            sprintf(filename,"vrlrpdfall.dat");
                            pdfarr(pG->time, -2.5 * maxvel, 2.5 * delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                            // Mass-weighted
                            sprintf(filename,"vrlrmasspdfall.dat");
                            pdfarr(pG->time, -2.5 * maxvel, 2.5 * delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                            // sprintf(filename,"vrlogpdfall.dat");
                            // pdfarr(pG->time, -1.0, 0.004, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 1);
                            // sprintf(filename,"vrlogmasspdfall.dat");
                            // pdfarr(pG->time, -1.0, 0.004, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 1);
                            // Velocity (radial) as a function of density
                            // PDF of radial velocity as a function of cell density
                            sprintf(filename,"vrdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 1, 1);
                            // PDF of radial velocity by LOS circumcluster density
                            sprintf(filename,"vrsigmadenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 2, 1);
                        
                            // PDFs of outflowing gas -  Density and velocity distributions but only for gas that is strictly outflowing
                            /* for (i = 1; i <= gnx1*gnx2*gnx3; i++) {
                                if (gblpr[i-1] < 0) {
                                    gblpr[i-1] = 0.0;
                                    gblden[i-1] = 0.0;
                                }
                            }
                            // Outflowing Velocity (radial) as a function of density
                            sprintf(filename,"ofdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, gnx1*gnx2*gnx3, 0, 1);
                            sprintf(filename,"vrofdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 1, 1);
                            sprintf(filename,"vrofsigmadenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 2, 1); */
                        }
                        
                        // Same but for infalling
                        /* expr = "Vr";
                        outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                        expr = "d";
                        outprocarr(pD, pG, allden, gblden, xcomvec, expr);
                        if (myID_Comm_world == 0) {
                            // PDFs of infalling gas
                            for (i = 1; i <= gnx1*gnx2*gnx3; i++) {
                                if (gblpr[i-1] >= 0) {
                                    gblpr[i-1] = 0.0;
                                    gblden[i-1] = 0.0;
                                }
                            }
                            // Infalling Velocity (radial) as a function of density
                            sprintf(filename,"ifdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, gnx1*gnx2*gnx3, 0, 1);
                            sprintf(filename,"vrifdenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 1, 1);
                            sprintf(filename,"vrifsigmadenpdfall.dat");
                            pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblpr, gnx1*gnx2*gnx3, 2, 1);
                        } */
                        
                        // PDF of x-velocity
                        /* expr = "V1";
                         outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                         if (myID_Comm_world == 0) {
                         sprintf(filename,"vxpdfall.dat");
                         pdfarr(pG->time, -1.0 * maxvel, delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 0, 0);
                         sprintf(filename,"vxmasspdfall.dat");
                         pdfarr(pG->time, -1.0 * maxvel, delvelr, npdf, filename, gblpr, gblden, gnx1*gnx2*gnx3, 1, 0);
                         } */
                        
                    }
                }
            }
            
            free_1d_array(allden);
            free_1d_array(gblden);
            free_1d_array(allpr);
            free_1d_array(gblpr);
            free_1d_array(gblsurfouterpdf);
            free_1d_array(gblanglepdf);
            free_1d_array(gblvsurfpdf);
            
        }
        
        // Flag indicates whether to output flux of several quantities described below at radial shells away from either the box center or the
        // stellar center of mass. 0 for no output or else fluxrad_out defines the number of radial bins
        if (fluxrad_out > 0) {
            
            // Calculate the density pdf in the real grid
            minden = log10(rho_floor);
            maxden = log10(rho_cloud) + 6.0;
            delden = (maxden - minden) / ((double) (npdf) - 1.0);
            
            // Allocate arrays as a function of theta, phi and shell radius for different flux quantities
            if ((gblxflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblyflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblzflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblxfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblyfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblzfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gbldflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((gblarea=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclxflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclyflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclzflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclxfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclyfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclzfluxpos=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lcldflux=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((lclarea=(Real*)calloc_1d_array(ntheta*nphi*nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for z-velocity grid\n");
            }
            if ((allrfluxes=(Real*)calloc_1d_array(nrs*ncoms*nfluxes,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for all fluxes\n");
            }
            if ((rflux=(Real*)calloc_1d_array(nrs,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for all fluxes\n");
            }
            if ((gblden=(Real*)calloc_1d_array(ntheta*nphi,sizeof(Real)))==NULL) {
                ath_error("[radpargrav]: Error allocating memory for PDF grid\n");
            }
            
            for (nl=0; nl<pM->NLevels; nl++) {
                for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++) {
                    if (pM->Domain[nl][nd].Grid != NULL) {
                        
                        for (nc = 1; nc <= ncoms; nc++) {
                            
                            // Radii at which to evaluate the flux in spheres and the origin defined by xcomvec
                            // nc = 1 is box center, and nc = 0 is stellar center of mass
                            if (nc == 1) {
                                rfluxmin = 0.01 * rcloud;
                                rfluxmax = 2.0 * rcloud - pG->dx1;
                                xcomvec[0] = 0.0;
                                xcomvec[1] = 0.0;
                                xcomvec[2] = 0.0;
                            } else {
                                rfluxmin = 0.01 * rcloud;
                                rfluxmax = 2.0 * rcloud - pG->dx1 - rcoms;
                                xcomvec[0] = xcoms;
                                xcomvec[1] = ycoms;
                                xcomvec[2] = zcoms;
                            }
                            delrflux = (rfluxmax - rfluxmin) / (nrs - 1.0);
                            krc = (int)round(1 + (rcloud - rfluxmin) / delrflux);
                            
                            // Initialization of vector containing output fluxes as a function of radius
                            for (k=1;k<=nrs;k++) {
                                rflux[k-1] = rfluxmin + (k-1) * delrflux;
                                for (nf = 1; nf <= nfluxes; nf++) {
                                    fluxind = (nc - 1) * nfluxes * nrs + (k - 1) * nfluxes + nf - 1;
                                    allrfluxes[fluxind] = 0.0;
                                    // ath_pout(0, "fluxind = %d of %d\n", fluxind, nrs * nfluxes * ncoms);
                                }
                                fluxind = (nc - 1) * nfluxes * nrs + (k - 1) * nfluxes;
                                allrfluxes[fluxind] = rflux[k - 1];
                            }
                            
                            // ath_pout(0, "Initializing density\n");
                            
                            // X Momentum Flux
                            expr = "d";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            // Evaluation of the density field in spheres around the box centre
                            loclcircfield(pG, rflux, lcldflux, lclarea, xcomvec, ntheta, nrs, 0, expr);
                            // ath_pout(0, "Density initialized\n");
                            expr = "M1";
                            // Evaluation of the px field in spheres around the box centre
                            loclcircfield(pG, rflux, lclxflux, lclarea, xcomvec, ntheta, nrs, 0, expr);
                            loclcircfield(pG, rflux, lclxfluxpos, lclarea, xcomvec, ntheta, nrs, 1, expr);
                            // ath_pout(0, "Momentum-1 initialized\n");
                            // Y Momentum Flux
                            expr = "M2";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            loclcircfield(pG, rflux, lclyflux, lclarea, xcomvec, ntheta, nrs, 0, expr);
                            loclcircfield(pG, rflux, lclyfluxpos, lclarea, xcomvec, ntheta, nrs, 2, expr);
                            // ath_pout(0, "Momentum-2 initialized\n");
                            // Z Momentum Flux
                            expr = "M3";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            loclcircfield(pG, rflux, lclzflux, lclarea, xcomvec, ntheta, nrs, 0, expr);
                            loclcircfield(pG, rflux, lclzfluxpos, lclarea, xcomvec, ntheta, nrs, 3, expr);
                            // ath_pout(0, "Momentum-3 initialized\n");
                            
                            
#ifdef MPI_PARALLEL
                            // Broadcast all arrays
                            mpierr = MPI_Reduce(lcldflux, gbldflux, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            // mpierr = MPI_Reduce(lclarea, gblarea, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclxflux, gblxflux, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclyflux, gblyflux, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclzflux, gblzflux, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclxfluxpos, gblxfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclyfluxpos, gblyfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclzfluxpos, gblzfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif
                            
                            if (myID_Comm_world == 0) {
                                
                                // Save fluxes in this order: Radius, Radiative flux, Tau, Mass Flux, Rho, Fr, Rho * Fr, Px, Py, Pz, Pr, Prout, Phi, dPhidr, RhodPhidr
                                
                                for (k = 1; k <= nrs; k++) {
                                    for (i = 1; i <= ntheta; i++) {
                                        for (j = 1; j <= nphi; j++) {
                                            iden = (k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                            tempflux = (gblxfluxpos[iden] + gblyfluxpos[iden] + gblzfluxpos[iden]) / rflux[k - 1];
                                            temprflux = tempflux * tempflux;
                                            fluxind = (nc - 1) * nfluxes * nrs + (k - 1) * nfluxes;
                                            allrfluxes[fluxind + 4] = allrfluxes[fluxind + 4] + lclarea[iden] * gbldflux[iden] / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // Rho
                                            allrfluxes[fluxind + 3] = allrfluxes[fluxind + 3] + lclarea[iden] * tempflux; // Mass Flux
                                            tempthetaflux = 0.0;
                                            for (kk = 1; kk <= k; kk++) {
                                                tempthetaflux = tempthetaflux + gbldflux[(kk - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] * delrflux;
                                            }
                                            if (k == nrs) {
                                                // ath_pout(0, "i = %d, j = %d, rhokappadr = %e\n", i, j, tempthetaflux * kappa_IR);
                                            }
                                            allrfluxes[fluxind + 2] = allrfluxes[fluxind + 2] + exp(-1.0 * tempthetaflux * kappa_IR) * lclarea[iden] / (4.0 * rflux[k - 1] * rflux[k - 1] * PI);
                                            allrfluxes[fluxind + 7] = allrfluxes[fluxind + 7] + (gblxflux[iden] * gblxfluxpos[iden] + gblxflux[iden] * gblyfluxpos[iden] + gblxflux[iden] * gblzfluxpos[iden]) * lclarea[iden] / (rflux[k - 1] * gbldflux[iden]); // Px
                                            allrfluxes[fluxind + 8] = allrfluxes[fluxind + 8] + (gblyflux[iden] * gblxfluxpos[iden] + gblyflux[iden] * gblyfluxpos[iden] + gblyflux[iden] * gblzfluxpos[iden]) * lclarea[iden] / (rflux[k - 1] * gbldflux[iden]); // Py
                                            allrfluxes[fluxind + 9] = allrfluxes[fluxind + 9] + (gblzflux[iden] * gblxfluxpos[iden] + gblzflux[iden] * gblyfluxpos[iden] + gblzflux[iden] * gblzfluxpos[iden]) * lclarea[iden] / (rflux[k - 1] * gbldflux[iden]); // Pz
                                            allrfluxes[fluxind + 10] = allrfluxes[fluxind + 10] + temprflux * lclarea[iden] / gbldflux[iden]; // Pr
                                            if (tempflux > 0) {
                                                allrfluxes[fluxind + 11] = allrfluxes[fluxind + 11] + temprflux * lclarea[iden] / gbldflux[iden]; // Prout
                                            }
                                        }
                                    }
                                }
                                // ath_pout(1, "Found all fluxes\n");
                                
                            }
                            // ath_error("Break\n");

                            
#ifdef RADIATION
                            // X Radiation Flux
                            expr = "F1";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            loclcircfield(pG, rflux, lclxfluxpos, lclarea, xcomvec, ntheta, nrs, 1, expr);
                            // ath_pout(0, "Flux-1 initialized\n");
                            
                            // Y Radiation Flux
                            expr = "F2";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            loclcircfield(pG, rflux, lclyfluxpos, lclarea, xcomvec, ntheta, nrs, 2, expr);
                            // ath_pout(0, "Flux-2 initialized\n");
                            
                            // Z Radiation Flux
                            expr = "F3";
                            // outprocarr(pD, pG, allpr, gblpr, xcomvec, expr);
                            loclcircfield(pG, rflux, lclzfluxpos, lclarea, xcomvec, ntheta, nrs, 3, expr);
                            // ath_pout(0, "Flux-3 initialized\n");
#endif
                            
                            // Gravitational field
                            expr = "Phi";
                            loclcircfield(pG, rflux, lclxflux, lclarea, xcomvec, ntheta, nrs, 0, expr);
                            
                            
#ifdef MPI_PARALLEL
                            mpierr = MPI_Reduce(lclxfluxpos, gblxfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclyfluxpos, gblyfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclzfluxpos, gblzfluxpos, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            mpierr = MPI_Reduce(lclxflux, gblxflux, nphi * ntheta * nrs, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
                            // ath_pout(0, "All arrays broadcast\n");
#endif
                            
                            // Fluxes related to the gradient of the potential
                            if (myID_Comm_world == 0) {
                                
                                for (k = 1; k <= nrs; k++) {
                                    for (i = 1; i <= ntheta; i++) {
                                        for (j = 1; j <= nphi; j++) {
                                            iden = (k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                            idenp = k * nphi * ntheta + (i - 1) * nphi + j - 1;
                                            idenm = (k - 2) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                            fluxind = (nc - 1) * nfluxes * nrs + (k - 1) * nfluxes;
                                            tempflux = (gblxfluxpos[iden] + gblyfluxpos[iden] + gblzfluxpos[iden]) / rflux[k - 1];
#ifdef RADIATION
                                            allrfluxes[fluxind + 1] = allrfluxes[fluxind + 1] + lclarea[iden] * tempflux; // Radiative Flux
                                            allrfluxes[fluxind + 5] = allrfluxes[fluxind + 5] + lclarea[iden] * tempflux / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // Fr
                                            allrfluxes[fluxind + 6] = allrfluxes[fluxind + 6] + lclarea[iden] * tempflux * gbldflux[iden] / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // Rho Fr
                                            if (k == 1) {
                                                temprflux = (gblxflux[idenp] - gblxflux[iden]) / delrflux;
                                            } else if (k == nrs) {
                                                temprflux = (gblxflux[iden] - gblxflux[idenm]) / delrflux;
                                            } else {
                                                temprflux = (gblxflux[idenp] - gblxflux[idenm]) / (2.0 * delrflux);
                                            }
#endif
                                            allrfluxes[fluxind + 14] = allrfluxes[fluxind + 14] + lclarea[iden] * gbldflux[iden] * temprflux / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // Rho DPhi/Dr
                                            allrfluxes[fluxind + 13] = allrfluxes[fluxind + 13] + lclarea[iden] * temprflux / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // DPhi/Dr
                                            allrfluxes[fluxind + 12] = allrfluxes[fluxind + 12] + lclarea[iden] * gblxflux[iden] / (4.0 * rflux[k - 1] * rflux[k - 1] * PI); // Phi
                                            // Calculation of density on shells, where the outward radial velocity is greater than zero
                                            if (tempflux > 0) {
                                                gblyflux[iden] = gbldflux[iden];
                                            } else {
                                                gblyflux[iden] = 0.0;
                                            }
                                        }
                                    }
                                }
                                
                                // PDFs of gas density for outflowing gas only
                                // ath_pout(1, "Found all fluxes\n");
                                // Density pdf of outflowing gas on shell at maximum radius
                                /* for (i = 1; i <= ntheta; i++) {
                                    for (j = 1; j <= nphi; j++) {
                                        idenp = (nrs - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                        iden = (i - 1) * nphi + j - 1;
                                        gblden[iden] = gblyflux[idenp];
                                    }
                                }
                                sprintf(filename,"denpdfof_rmax_r%d.dat", nc - 1);
                                pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, ntheta*nphi, 0, 1);
                                
                                // Density pdf of all gas on shell at maximum radius
                                for (i = 1; i <= ntheta; i++) {
                                    for (j = 1; j <= nphi; j++) {
                                        idenp = (nrs - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                        iden = (i - 1) * nphi + j - 1;
                                        gblden[iden] = gbldflux[idenp];
                                    }
                                }
                                sprintf(filename,"denpdfall_rmax_r%d.dat", nc - 1);
                                pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, ntheta*nphi, 0, 1);
                                
                                // Density pdf of outflowing gas on shell at cloud radius
                                for (i = 1; i <= ntheta; i++) {
                                    for (j = 1; j <= nphi; j++) {
                                        idenp = (krc - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                        iden = (i - 1) * nphi + j - 1;
                                        gblden[iden] = gblyflux[idenp];
                                    }
                                }
                                sprintf(filename,"denpdfof_r0_r%d.dat", nc - 1);
                                pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, ntheta*nphi, 0, 1);
                                
                                // Density pdf of all gas on shell at cloud radius
                                for (i = 1; i <= ntheta; i++) {
                                    for (j = 1; j <= nphi; j++) {
                                        idenp = (krc - 1) * nphi * ntheta + (i - 1) * nphi + j - 1;
                                        iden = (i - 1) * nphi + j - 1;
                                        gblden[iden] = gbldflux[idenp];
                                    }
                                }
                                sprintf(filename,"denpdfall_r0_r%d.dat", nc - 1);
                                pdfarr(pG->time, minden, delden, npdf, filename, gblden, gblden, ntheta*nphi, 0, 1); */
                            }
                            
                        }
                        
                        if (myID_Comm_world == 0) {
                            
                            // ath_pout(0, "Beginning first dump of stars and fluxes\n");
                            
                            // Output relevant flux quantities and stellar centres of mass
                            mode = (pG->time == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                            sprintf(filename,"../com_star.dat");
#else
                            sprintf(filename,"com_star.dat");
#endif
                            if((hp = fopen(filename,mode)) == NULL) {
                                ath_error("[radpargrav]: Unable to open density dump file\n");
                                return;
                            }
                            
                            fprintf(hp,"\n%lf %e %e %e %e\n", pG->time, xcoms, ycoms, zcoms, rcoms);
                            fclose(hp);
                            
                            for (nc = 1; nc <= ncoms; nc++) {
                                
                                mode = (pG->time == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
                                sprintf(filename,"../fluxes_r%d.dat", nc - 1);
#else
                                sprintf(filename,"fluxes_r%d.dat", nc - 1);
#endif
                                if((hp = fopen(filename,mode)) == NULL) {
                                    ath_error("[radpargrav]: Unable to open density dump file\n");
                                    return;
                                }
                                
                                // Output fluxes in this order: Radius, Radiative flux, Tau, Mass Flux, Rho, Fr, Rho * Fr, Px, Py, Pz, Pr, Prout, Phi, dPhidr, RhodPhidr
                                for (k = 1; k <= nrs; k++) {
                                    
                                    iden = (nc - 1) * nfluxes * nrs + (k - 1) * nfluxes;
                                    
                                    fprintf(hp,"\n%lf %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", pG->time, allrfluxes[iden], allrfluxes[iden + 1], allrfluxes[iden + 2], allrfluxes[iden + 3], allrfluxes[iden + 4], allrfluxes[iden + 5], allrfluxes[iden + 6], allrfluxes[iden + 7], allrfluxes[iden + 8], allrfluxes[iden + 9], allrfluxes[iden + 10], allrfluxes[iden + 11], allrfluxes[iden + 12], allrfluxes[iden + 13], allrfluxes[iden + 14]);
                                    
                                }
                                
                                fclose(hp);
                                
                            }
                            
                            // ath_pout(0, "Ending first dump of stars and fluxes\n");
                            
                        }
                        
                        // ath_error("Break\n");
                        
                    }
                }
            }
            free_1d_array(gblxflux);
            free_1d_array(gblyflux);
            free_1d_array(gblzflux);
            free_1d_array(gblxfluxpos);
            free_1d_array(gblyfluxpos);
            free_1d_array(gblzfluxpos);
            free_1d_array(gbldflux);
            free_1d_array(gblarea);
            free_1d_array(lclxflux);
            free_1d_array(lclyflux);
            free_1d_array(lclzflux);
            free_1d_array(lclxfluxpos);
            free_1d_array(lclyfluxpos);
            free_1d_array(lclzfluxpos);
            free_1d_array(lcldflux);
            free_1d_array(lclarea);
            free_1d_array(allrfluxes);
            free_1d_array(rflux);
            free_1d_array(gblden);
            
        }
        
        den_counter++;
        den_timec = den_timec + rho_fdt;
            
    }
    
    if (isnan(pG->dt)) ath_error("[radpargrav]:  Time step is NaN!");
    
    if (idrive == 0) {  /* driven turbulence */
        /* Integration has already been done, but time not yet updated */
        newtime = pG->time + pG->dt;
        
#ifndef IMPULSIVE_DRIVING
        /* Drive on every time step */
        perturb(pG, pG->dt);
#endif /* IMPULSIVE_DRIVING */
        
        if (newtime >= (tdrive+dtdrive)) {
            /* If we start with large time steps so that tdrive would get way
             * behind newtime, this makes sure we don't keep generating after
             * dropping down to smaller time steps */
            while ((tdrive+dtdrive) <= newtime) tdrive += dtdrive;
            
#ifdef IMPULSIVE_DRIVING
            /* Only drive at intervals of dtdrive */
            perturb(pG, dtdrive);
#endif /* IMPULSIVE_DRIVING */
            
            /* Compute new spectrum after dtdrive.  Putting this after perturb()
             * means we won't be applying perturbations from a new power spectrum
             * just before writing outputs.  At the very beginning, we'll go a
             * little longer before regenerating, but the energy injection rate
             * was off on the very first timestep anyway.  When studying driven
             * turbulence, all we care about is the saturated state. */
            generate();
        }
    }
    
    return;
}

void Userwork_in_rad_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  Userwork_in_loop(pM);
  
  return;
}

/*==============================================================================
 * PRIVATE FUNCTIONS
 *============================================================================*/

/*------------------------------------------------------------------------------
 *  Function hst_*
 *
 *  Dumps to history file
 *  "History" variables are usually volume integrals: e.g., 
 *    S = \Sum_{ijk} q[k][j][i]*dx1[i]*dx2[j]*dx3[k].  Note that history 
 *    variables are TOTAL, VOLUME INTEGRATED quantities, NOT volume averages 
 *    (just divide by the Domain volume to compute the latter).  With MPI, the
 *    sum is performed over all Grids in the Domain.
 */

/* Dump kinetic energy in perturbations */
static Real hst_dEk(const GridS *pG,const int i,const int j,const int k)
{
  /* The kinetic energy in perturbations is 0.5*d*V^2 */
  return 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 +
              pG->U[k][j][i].M2*pG->U[k][j][i].M2 +
              pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
}

/* Dump magnetic energy in perturbations */
static Real hst_dEb(const GridS *pG,const int i,const int j,const int k)
{
  /* The magnetic energy in perturbations is 0.5*B^2 - 0.5*B0^2 */
#ifdef MHD
  return 0.5*((pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
               pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
               pG->U[k][j][i].B3c*pG->U[k][j][i].B3c)-B0*B0);
#else /* MHD */
  return 0.0;
#endif /* MHD */
}

// Magnetic energy in perturbations within the original cloud radius
static Real hst_dEb_one(const GridS *pG,const int i,const int j,const int k)
{
    Real x1,x2,x3,tmp;
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    tmp = (SQR(x1) + SQR(x2) + SQR(x3));
    /* The magnetic energy in perturbations is 0.5*B^2 - 0.5*B0^2 */
#ifdef MHD
    if (tmp<=1.0*rcloud*rcloud) {
        return 0.5*((pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
                     pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
                     pG->U[k][j][i].B3c*pG->U[k][j][i].B3c)-B0*B0);
    } else {
        return 0.0;
    }
#else /* MHD */
    return 0.0;
#endif /* MHD */
}

/* Dump mass in star particles */
static Real hst_mass_stars(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real x1,x2,x3;
  
  pList = pG->Lstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    /* If particle in this cell */
    if (i==ip && j==jp && k==kp) {
      return pStar->m/(pG->dx1*pG->dx2*pG->dx3);
    }
    pList = pList->next;
  }
  return 0.0;
}

/* Dump mass in gas */
static Real hst_mass_gas(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
    
  /* The gravitational energy is d*Phi, except where particles are, so check to
   * see if there is a particle. */
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1)
      return 0.0;
    pList = pList->next;
  }
  return pG->U[k][j][i].d;
}

/* Dump gravitational potential energy in gas plus stars */
static Real hst_Egrav_tot(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real tmp,x1,x2,x3,W1[3],W2[3],W3[3];
  
  /* The gravitational energy is d*Phi, except where particles are, so check to
   * see if there is a particle. */
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    /* If within 1 zone of a star particle */
    if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      /* Compute TSC weights   */
      tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
      W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
      tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
      W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
      tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
      W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
      
      /* Return Phi times weighted particle density */
      tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
      return tmp*pG->Phi[k][j][i]*W1[i-(ip-1)]*W2[j-(jp-1)]*W3[k-(kp-1)];
    }
    pList = pList->next;
  }
  return pG->U[k][j][i].d*pG->Phi[k][j][i];
}

/* Dump gravitational potential energy in gas */
static Real hst_Egrav_gas(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  
  /* The gravitational energy is d*Phi, except where particles are, so check to
   * see if there is a particle. */
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1)
      return 0.0;
    pList = pList->next;
  }
  return pG->U[k][j][i].d*pG->Phi[k][j][i];
}

/* Dump gravitational potential energy in stars */
static Real hst_Egrav_stars(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real tmp,x1,x2,x3,W1[3],W2[3],W3[3];
  
  /* The gravitational energy is d*Phi, except where particles are, so check to
   * see if there is a particle. */
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    /* If within 1 zone of a star particle */
    if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      /* Compute TSC weights   */
      tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
      W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
      tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
      W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
      tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
      W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
      
      /* Return Phi times weighted particle density */
      tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
      return tmp*pG->Phi[k][j][i]*W1[i-(ip-1)]*W2[j-(jp-1)]*W3[k-(kp-1)];
    }
    pList = pList->next;
  }
  return 0.0;
}

/* Dump rate of work done by gravitational field on gas */
static Real hst_Wgrav_gas(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real x1,x2,x3,f1,f2,f3;
  
  /* The gravitational work is -d*(v dot Phi), except where particles are, so 
   * first check to see if there is a particle. */
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    if (abs(i-ip)<=1 && abs(j-jp)<=1 && abs(k-kp)<=1)
      return 0.0;
    pList = pList->next;
  }
  /* Otherwise, calculate force using centered-difference approximation */
  f1 = -0.5*(pG->Phi[k][j][i+1]-pG->Phi[k][j][i-1])/pG->dx1;
  f2 = -0.5*(pG->Phi[k][j+1][i]-pG->Phi[k][j-1][i])/pG->dx2;
  f3 = -0.5*(pG->Phi[k+1][j][i]-pG->Phi[k-1][j][i])/pG->dx3;
  return (pG->U[k][j][i].M1*f1 +
          pG->U[k][j][i].M2*f2 +
          pG->U[k][j][i].M3*f3)/pG->U[k][j][i].d;
}

/* Dump rate of work done by gravitational field on stars */
static Real hst_Wgrav_stars(const GridS *pG,const int i,const int j,const int k)
{
  StarParListS *pList=pG->Gstars;
  StarParS *pStar=NULL;
  int ip,jp,kp;
  Real tmp,x1,x2,x3,f1,f2,f3,W1[3],W2[3],W3[3];
  
  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
    cc_ijk(pG,pStar->x1,pStar->x2,pStar->x3,&ip,&jp,&kp);
    /* If particle in this cell */
    if (i==ip && j==jp && k==kp) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      /* Compute TSC weights   */
      tmp = 0.5 + (pStar->x1 - x1)/pG->dx1;
      W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2]-tmp+0.5;  W1[1] = 1.0-W1[0]-W1[2];
      tmp = 0.5 + (pStar->x2 - x2)/pG->dx2;
      W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2]-tmp+0.5;  W2[1] = 1.0-W2[0]-W2[2];
      tmp = 0.5 + (pStar->x3 - x3)/pG->dx3;
      W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2]-tmp+0.5;  W3[1] = 1.0-W3[0]-W3[2];
      
      /* Return v dot grad(Phi) times weighted particle density */
      tmp = pStar->m/(pG->dx1*pG->dx2*pG->dx3);
      f1 = -0.5*((pG->Phi[kp][jp][ip  ] - pG->Phi[kp][jp][ip-2])*W1[0] +
                 (pG->Phi[kp][jp][ip+1] - pG->Phi[kp][jp][ip-1])*W1[1] +
                 (pG->Phi[kp][jp][ip+2] - pG->Phi[kp][jp][ip  ])*W1[2])/pG->dx1;
      f2 = -0.5*((pG->Phi[kp][jp  ][ip] - pG->Phi[kp][jp-2][ip])*W2[0] +
                 (pG->Phi[kp][jp+1][ip] - pG->Phi[kp][jp-1][ip])*W2[1] +
                 (pG->Phi[kp][jp+2][ip] - pG->Phi[kp][jp  ][ip])*W2[2])/pG->dx2;
      f3 = -0.5*((pG->Phi[kp  ][jp][ip] - pG->Phi[kp-2][jp][ip])*W3[0] +
                 (pG->Phi[kp+1][jp][ip] - pG->Phi[kp-1][jp][ip])*W3[1] +
                 (pG->Phi[kp+2][jp][ip] - pG->Phi[kp  ][jp][ip])*W3[2])/pG->dx3;
      
      return tmp*(f1*pStar->v1 + f2*pStar->v2 + f3*pStar->v3);
    }
    pList = pList->next;
  }
  return 0.0;
}

/* Dump effective squared radius weighted by the density */
static Real hst_Reffone_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = (SQR(x1) + SQR(x2) + SQR(x3));
  if (tmp<=1.0*rcloud*rcloud) {
    return tmp*pG->U[k][j][i].d;
  }
  return 0.0;
}

/* Dump effective squared radius weighted by the density within original cloud radius */
static Real hst_rhoeffone_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = (SQR(x1) + SQR(x2) + SQR(x3));
  if (tmp<=1.0*rcloud*rcloud) {
    return pG->U[k][j][i].d;
  }
  return 0.0;
}

/* Dump density-weighted density within original cloud radius */
static Real hst_Meffone_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = (SQR(x1) + SQR(x2) + SQR(x3));
  if (tmp<=1.0*rcloud*rcloud) {
    return pG->U[k][j][i].d*pG->U[k][j][i].d;
  }
  return 0.0;
}

/* Dump density-weighted r^2*/
static Real hst_Reff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = (SQR(x1) + SQR(x2) + SQR(x3));
  return tmp*pG->U[k][j][i].d;
}

// Effective positions used to calculated principal directions of gas distribution
/* Dump effective x position weighted by the density */
static Real hst_xeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x1*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective y position weighted by the density */
static Real hst_yeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x2*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective z position weighted by the density */
static Real hst_zeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x3*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective x^2 position weighted by the density */
static Real hst_xxeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x1*x1*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective y^2 position weighted by the density */
static Real hst_yyeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x2*x2*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective z^2 position weighted by the density */
static Real hst_zzeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x3*x3*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective x*y position weighted by the density */
static Real hst_xyeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x1*x2*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective x*z position weighted by the density */
static Real hst_xzeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x1*x3*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective y*z position weighted by the density */
static Real hst_yzeff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp;
    
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = x2*x3*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump effective mass weighted by the density */
static Real hst_Meff_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp;
    
  tmp = pG->U[k][j][i].d*pG->U[k][j][i].d;
    
  return tmp;
}

/* Dump radial component of the outward momentum flux */
static Real hst_Mr_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp,rmf=0.0;
  
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = (x1*pG->U[k][j][i].M1 +
         x2*pG->U[k][j][i].M2 +
         x3*pG->U[k][j][i].M3)/sqrt(SQR(x1) + SQR(x2) + SQR(x3));

  /* Inner x1-boundary */
  if (i==pG->is && pG->lx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
  
  /* Outer x1-boundary */
  if (i==pG->ie && pG->rx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
  
  /* Inner x2-boundary */
  if (j==pG->js && pG->lx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
  
  /* Outer x2-boundary */
  if (j==pG->je && pG->rx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
  
  /* Inner x3-boundary */
  if (k==pG->ks && pG->lx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  
  /* Outer x3-boundary */
  if (k==pG->ke && pG->rx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  
  return rmf;
}

/* Dump outward free energy flux */
/* NOTE:  This is actually an energy density flux, since history dumps are 
 *   volume integrals */
static Real hst_Efree_out(const GridS *pG,const int i,const int j,const int k)
{
  Real x1,x2,x3,tmp,ef=0.0;

  /* Free Energy = Kinetic + Thermal + Gravitational */
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  tmp = 0.5*(SQR(pG->U[k][j][i].M1) +
             SQR(pG->U[k][j][i].M2) +
             SQR(pG->U[k][j][i].M3))/SQR(pG->U[k][j][i].d) +
    1.5*Iso_csound2 + pG->Phi[k][j][i];
  
  /* Inner x1-boundary */
  if (i==pG->is && pG->lx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    ef -= tmp*pG->U[k][j][i].M1*pG->dt/pG->dx1;
  
  /* Outer x1-boundary */
  if (i==pG->ie && pG->rx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    ef += tmp*pG->U[k][j][i].M1*pG->dt/pG->dx1;
  
  /* Inner x2-boundary */
  if (j==pG->js && pG->lx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    ef -= tmp*pG->U[k][j][i].M2*pG->dt/pG->dx2;
  
  /* Outer x2-boundary */
  if (j==pG->je && pG->rx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    ef += tmp*pG->U[k][j][i].M2*pG->dt/pG->dx2;
  
  /* Inner x3-boundary */
  if (k==pG->ks && pG->lx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    ef -= tmp*pG->U[k][j][i].M3*pG->dt/pG->dx3;
  
  /* Outer x3-boundary */
  if (k==pG->ke && pG->rx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    ef += tmp*pG->U[k][j][i].M3*pG->dt/pG->dx3;
  
  return ef;
}

/* Dump outward mass flux */
static Real hst_Mass_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].d);
    
  /* Inner x1-boundary */
  if (i==pG->is && pG->lx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
    
  /* Outer x1-boundary */
  if (i==pG->ie && pG->rx1_id==-1 &&
      j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
    rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
    
  /* Inner x2-boundary */
  if (j==pG->js && pG->lx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
    
  /* Outer x2-boundary */
  if (j==pG->je && pG->rx2_id==-1 &&
      i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
    rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
    
  /* Inner x3-boundary */
  if (k==pG->ks && pG->lx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
    
  /* Outer x3-boundary */
  if (k==pG->ke && pG->rx3_id==-1 &&
      i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
    rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
    
  return rmf;
    
}

// If we define scalars to track mass in original locations, then calculate mass flux outwards from these different passive scalars
#if (NSCALARS > 0)
static Real hst_MassS0_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].s[0]);
    
  /* Inner x1-boundary */
  if (tmp > 0.0) {
    if (i==pG->is && pG->lx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
    
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
    
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
    
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
    
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
    
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  }
    
  return rmf;
    
}

static Real hst_MassS1_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].s[1]);
    
  /* Inner x1-boundary */
  if (tmp > 0.0) {
    if (i==pG->is && pG->lx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
        
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  }
    
  return rmf;
    
}

static Real hst_MassS2_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].s[2]);
    
  /* Inner x1-boundary */
  if (tmp > 0.0) {
    if (i==pG->is && pG->lx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
        
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  }
    
  return rmf;
    
}

static Real hst_MassS3_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].s[3]);
    
  /* Inner x1-boundary */
  if (tmp > 0.0) {
    if (i==pG->is && pG->lx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
        
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  }
    
  return rmf;
    
}

static Real hst_MassS4_out(const GridS *pG,const int i,const int j,const int k)
{
  Real tmp,rmf=0.0;
    
  tmp = (pG->U[k][j][i].s[4]);
    
  /* Inner x1-boundary */
  if (tmp > 0.0) {
    if (i==pG->is && pG->lx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1 &&
	j>=pG->js && j<=pG->je && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M1/pG->U[k][j][i].d)*pG->dt/pG->dx1;
        
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf -= tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1 &&
	i>=pG->is && i<=pG->ie && k>=pG->ks && k<=pG->ke)
      rmf += tmp*(pG->U[k][j][i].M2/pG->U[k][j][i].d)*pG->dt/pG->dx2;
        
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf -= tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
        
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1 &&
	i>=pG->is && i<=pG->ie && j>=pG->js && j<=pG->je)
      rmf += tmp*(pG->U[k][j][i].M3/pG->U[k][j][i].d)*pG->dt/pG->dx3;
  }
    
  return rmf;
    
}
#endif

/* Dump outward radiation flux */
static Real hst_F_out(const GridS *pG,const int i,const int j,const int k)
{
  Real fout=0.0;

  /* Inner x1-boundary */
#ifdef RADIATION
  if (i==pG->is && pG->lx1_id==-1)
    fout -= pG->Urad[k][j][i].F1/pG->dx1;
  
  /* Outer x1-boundary */
  if (i==pG->ie && pG->rx1_id==-1)
    fout += pG->Urad[k][j][i].F1/pG->dx1;
  
  /* Inner x2-boundary */
  if (j==pG->js && pG->lx2_id==-1)
    fout -= pG->Urad[k][j][i].F2/pG->dx2;
  
  /* Outer x2-boundary */
  if (j==pG->je && pG->rx2_id==-1)
    fout += pG->Urad[k][j][i].F2/pG->dx2;
  
  /* Inner x3-boundary */
  if (k==pG->ks && pG->lx3_id==-1)
    fout -= pG->Urad[k][j][i].F3/pG->dx3;
 
  /* Outer x3-boundary */
  if (k==pG->ke && pG->rx3_id==-1)
    fout += pG->Urad[k][j][i].F3/pG->dx3;
#endif
  
  return fout;
}

/* Dump outward radiation flux correlated with density */
static Real hst_Frho_out(const GridS *pG,const int i,const int j,const int k)
{
    Real fout=0.0;
    
    /* Inner x1-boundary */
#ifdef RADIATION
    if (i==pG->is && pG->lx1_id==-1)
        fout -= pG->Urad[k][j][i].F1 * pG->U[k][j][i].d/pG->dx1;
    
    /* Outer x1-boundary */
    if (i==pG->ie && pG->rx1_id==-1)
        fout += pG->Urad[k][j][i].F1 * pG->U[k][j][i].d/pG->dx1;
    
    /* Inner x2-boundary */
    if (j==pG->js && pG->lx2_id==-1)
        fout -= pG->Urad[k][j][i].F2 * pG->U[k][j][i].d/pG->dx2;
    
    /* Outer x2-boundary */
    if (j==pG->je && pG->rx2_id==-1)
        fout += pG->Urad[k][j][i].F2 * pG->U[k][j][i].d/pG->dx2;
    
    /* Inner x3-boundary */
    if (k==pG->ks && pG->lx3_id==-1)
        fout -= pG->Urad[k][j][i].F3 * pG->U[k][j][i].d/pG->dx3;
    
    /* Outer x3-boundary */
    if (k==pG->ke && pG->rx3_id==-1)
        fout += pG->Urad[k][j][i].F3 * pG->U[k][j][i].d/pG->dx3;
#endif
    
    return fout;
}


/* Dump total source radiation */
static Real hst_jsrc(const GridS *pG,const int i,const int j,const int k)
{
  Real jsrc=0.0,rsq_src,tmp1,tmp2,tmp3,rsq;
  StarParListS *pList=NULL;
  StarParS *pStar=NULL;
  Real x1,x2,x3;
  
#ifdef RADIATION
  rsq_src = SQR(rsrc)/(2.0*log(2.0));
  tmp1 = pow(2.0*PI*rsq_src,-1.5);
  tmp2 = -0.5/rsq_src;

  pList = pG->Gstars;
  while (pList) {
    pStar = &(pList->starpar);
  
    tmp3 = tmp1*LuminosityFun(pStar->m,pStar->age);
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    rsq = SQR(x1-pStar->x1) + SQR(x2-pStar->x2) + SQR(x3-pStar->x3);
    jsrc += tmp3*exp(tmp2*rsq);
    pList = pList->next;
  }
#endif

  return jsrc;
}

/* Dump total radiation energy */
static Real hst_Er(const GridS *pG,const int i,const int j,const int k)
{
    Real Erout=0.0;
#ifdef RADIATION
  Erout = pG->Urad[k][j][i].Er;
#endif
    return Erout;
}

/* Dump total magnetic fields */
static Real hst_B3c(const GridS *pG,const int i,const int j,const int k)
{
    Real Bout=0.0;
#ifdef MHD
    Bout = pG->U[k][j][i].B3c;
#endif
    return Bout;
}


/*------------------------------------------------------------------------------
 *  Function usr_*
 *
 *  User-defined output expressions.  Surface densities are produced by using
 *  the ':' reduction operator in the input file, which produces a global 
 *  average along a given dimension.  For example, to obtain Sigma1, multiplying 
 *  by Lx and averaging is the same as multiplying the sum over all i of 
 *  rho(i,j,k) by dx=Lx/Nx.
 */

static Real log10d(const GridS *pG,const int i,const int j,const int k)
{
  return log10(pG->U[k][j][i].d);
}

static Real usr_Sigma1(const GridS *pG,const int i,const int j,const int k)
{
  return pG->U[k][j][i].d*Lx;
}

static Real usr_Sigma2(const GridS *pG,const int i,const int j,const int k)
{
  return pG->U[k][j][i].d*Ly;
}

static Real usr_Sigma3(const GridS *pG,const int i,const int j,const int k)
{
  return pG->U[k][j][i].d*Lz;
}


static Real const_absorption(const ConsS *pU)
{
  return kappa_IR;
}

static Real mass_dependent_luminosity(Real mass, Real age)
{
  return Psi*mass;
}

static Real rho(Real r)
{
  return (r <= rcloud) ? M_GMC/(FOUR_3RDS*PI*CUBE(rcloud)) : 0.0;
}

/*------------------------------------------------------------------------------
 *  Function diode_outflow_*
 *
 *  DIODE outflow boundary condition functions.
 */

static void diode_outflow_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][is-i] = pG->U[k][j][is];
        pG->U[k][j][is-i].M1 = MIN(pG->U[k][j][is-i].M1,0.0);
      }
    }
  }
}

static void diode_outflow_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pG->U[k][j][ie+i] = pG->U[k][j][ie];
        pG->U[k][j][ie+i].M1 = MAX(pG->U[k][j][ie+i].M1,0.0);
      }
    }
  }
}

static void diode_outflow_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][js-j][i] = pG->U[k][js][i];
        pG->U[k][js-j][i].M2 = MIN(pG->U[k][js-j][i].M2,0.0);
      }
    }
  }
}

static void diode_outflow_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[k][je+j][i] = pG->U[k][je][i];
        pG->U[k][je+j][i].M2 = MAX(pG->U[k][je+j][i].M2,0.0);
      }
    }
  }
}

static void diode_outflow_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ks-k][j][i] = pG->U[ks][j][i];
        pG->U[ks-k][j][i].M3 = MIN(pG->U[ks-k][j][i].M3,0.0);
      }
    }
  }
}

static void diode_outflow_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;
  
  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->U[ke+k][j][i] = pG->U[ke][j][i];
        pG->U[ke+k][j][i].M3 = MAX(pG->U[ke+k][j][i].M3,0.0);
      }
    }
  }
}


/*------------------------------------------------------------------------------
 *  Functions used to write arrays to file and calculate pdfs
 *
 *----------------------------------------------------------------------------*/

// Given a time, min-density, max-density, #of density bins in PDF, output file for PDF and integer array storing the PDF,
// outputs the PDF to file
static void intarrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, int *datarr)
{

  int i;
  char filename[24];
  char *mode;
  FILE *hp;

  mode = (tval == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
  sprintf(filename,"../");
  strcat(filename,filestr);
#else
  sprintf(filename,filestr);
#endif
  if((hp = fopen(filename,mode)) == NULL) {
    ath_error("[radpargrav]: Unable to open density dump file\n");
    return;
  }
  fprintf(hp,"\n%lf %lf %lf %d\n", tval, minval, delval, nvals);
  for (i=1; i<=nvals; i++) {
    fprintf(hp,"%d\t", datarr[i-1]);
  }
  fclose(hp);

}

// Same as intarrdatout but for floating point arrays - Probably a better way to do this
static void fltarrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *datarr)
{

  int i;
  char filename[24];
  char *mode;
  FILE *hp;

  mode = (tval == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
  sprintf(filename,"../");
  strcat(filename,filestr);
#else
  sprintf(filename,filestr);
#endif
  if((hp = fopen(filename,mode)) == NULL) {
    ath_error("[radpargrav]: Unable to open density dump file\n");
    return;
  }
  fprintf(hp,"\n%lf %lf %lf %d\n", tval, minval, delval, nvals);
  for (i=1; i<=nvals; i++) {
    fprintf(hp,"%e\t", datarr[i-1]);
  }
  fclose(hp);

}

// Same as fltarrdatout, but saves space by omitting values that are zero (below the lowest density
// and above the highest)
static void fltziparrdatout(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *datarr)
{

  int i, initind, maxind, initflag, maxflag;
  char filename[24];
  char *mode;
  FILE *hp;

  mode = (tval == 0.0) ? "w" : "a";
#ifdef MPI_PARALLEL
  sprintf(filename,"../");
  strcat(filename,filestr);
#else
  sprintf(filename,filestr);
#endif
  if((hp = fopen(filename,mode)) == NULL) {
    ath_error("[radpargrav]: Unable to open density dump file\n");
    return;
  }
  initind = 0;
  maxind = 0;
  maxflag = 0;
  initflag = 0;
  for (i = 1; i <=nvals; i++) {
    if ((datarr[i-1] != 0.0) && (initflag == 0)) {
      initind = i - 1;
      initflag = 1;
    }
    if ((datarr[i-1] != 0.0)) {
      maxind = i - 1;
    }
  }
  fprintf(hp,"\n%lf %lf %lf %d %d %d\n", tval, minval, delval, nvals, initind, maxind - initind + 1);
  for (i=initind; i<=maxind; i++) {
    fprintf(hp,"%e\t", datarr[i]);
  }
  fclose(hp);

}

// Given a time, min density, max density, # of bins in PDF, output file for PDF, real arrays storing the grid of density values, and the amount by which to
// augment the PDF for each value on that grid, outputs the PDF for densities on the grid. lflag for logarithmic bins, iflag for integer output PDF with
// max integer set by imax
static void pdfarr(Real tval, Real minval, Real delval, int nvals, char *filestr, Real *gridarr, Real *augarr, int imax, int iflag, int lflag)
{

  int i, iden;
  int *datarr=NULL;
  Real *rdatarr=NULL;

  if (iflag == 0) {
    if ((datarr=(int*)calloc_1d_array(nvals,sizeof(int)))==NULL) {
      ath_error("[radpargrav]: Error allocating memory for pdf\n");
    }
  } else {
    if ((rdatarr=(Real*)calloc_1d_array(nvals,sizeof(Real)))==NULL) {
      ath_error("[radpargrav]: Error allocating memory for pdf\n");
    }
  }

  for (i=1; i<=imax; i++) {
    if (lflag > 0) {
        if (gridarr[i-1] > 0.0) {
            iden = (int) ((log10(gridarr[i-1]) - minval) / delval);
        } else {
            iden = 0.0;
        }
    } else {
      iden = (int) ((gridarr[i-1] - minval) / delval);
    }
    if (iden >= nvals) {
      iden=nvals-1;
    }
    if (iden < 0) {
      iden=0;
    }
    if (iflag == 0) {
      datarr[iden]=datarr[iden]+1;
    } else if (iflag == 1) {
      rdatarr[iden]=rdatarr[iden]+augarr[i-1];
    } else {
        rdatarr[iden]=rdatarr[iden]+augarr[i-1]*augarr[i-1];
    }
  }

  if (iflag == 0) {
    intarrdatout(tval, minval, delval, nvals, filestr, datarr);
    free_1d_array(datarr);
  } else {
    fltziparrdatout(tval, minval, delval, nvals, filestr, rdatarr);
    free_1d_array(rdatarr);
  }

}

// Shortcut for transforming Athena grid values on pG to 1D array storing all grid values in gridarr. expr denotes
// which quantity to store in gridarr
static void outprocarr(const DomainS *pD, const GridS *pG, Real *locarr, Real *gridarr, Real *xcomvec, char *expr)
{

  int i, j, k, mpierr;
  Real tv, xd, yd, zd, xcoms, ycoms, zcoms, ddx, ddy, ddz;

  xcoms = xcomvec[0];
  ycoms = xcomvec[1];
  zcoms = xcomvec[2];

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	if (strcmp(expr,"d")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].d;
	} else if (strcmp(expr,"M1")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M1;
	} else if (strcmp(expr,"M2")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M2;
	} else if (strcmp(expr,"M3")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M3;
	} else if (strcmp(expr,"Mr")==0) {
	  cc_pos(pG, i, j, k, &xd, &yd, &zd);
	  tv = sqrt(SQR(xd - xcoms) + SQR(yd-ycoms) + SQR(zd-zcoms));
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = (pG->U[k][j][i].M1 * (xd-xcoms) + pG->U[k][j][i].M2 * (yd-ycoms) + pG->U[k][j][i].M3 * (zd-zcoms)) / tv;
	} else if (strcmp(expr,"V1")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
	} else if (strcmp(expr,"V2")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
	} else if (strcmp(expr,"V3")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
	} else if (strcmp(expr,"Vr")==0) {
	  cc_pos(pG, i, j, k, &xd, &yd, &zd);
	  tv = sqrt(SQR(xd - xcoms) + SQR(yd-ycoms) + SQR(zd-zcoms));
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = (pG->U[k][j][i].M1 * (xd-xcoms) + pG->U[k][j][i].M2 * (yd-ycoms) + pG->U[k][j][i].M3 * (zd-zcoms)) / (tv * pG->U[k][j][i].d);
	} else if (strcmp(expr,"Vt")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = sqrt(pG->U[k][j][i].M1*pG->U[k][j][i].M1+pG->U[k][j][i].M2*pG->U[k][j][i].M2+pG->U[k][j][i].M3*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
    } else if (strcmp(expr,"Mt")==0) {
        locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = sqrt(pG->U[k][j][i].M1*pG->U[k][j][i].M1+pG->U[k][j][i].M2*pG->U[k][j][i].M2+pG->U[k][j][i].M3*pG->U[k][j][i].M3);
#ifdef RADIATION
    } else if (strcmp(expr,"F1")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->Urad[k][j][i].F1;
	} else if (strcmp(expr,"F2")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->Urad[k][j][i].F2;
	} else if (strcmp(expr,"F3")==0) {
	  locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->Urad[k][j][i].F3;
    } else if (strcmp(expr,"Fr")==0) {
        cc_pos(pG, i, j, k, &xd, &yd, &zd);
        tv = sqrt(SQR(xd - xcoms) + SQR(yd-ycoms) + SQR(zd-zcoms));
        locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = (pG->Urad[k][j][i].F1 * (xd-xcoms) + pG->Urad[k][j][i].F2 * (yd-ycoms) + pG->Urad[k][j][i].F3 * (zd-zcoms)) / tv;
#endif
    } else if (strcmp(expr,"Phi")==0) {
        locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = pG->Phi[k][j][i];
    } else if (strcmp(expr,"dPhidr")==0) {
        cc_pos(pG, i, j, k, &xd, &yd, &zd);
        tv = sqrt(SQR(xd - xcoms) + SQR(yd-ycoms) + SQR(zd-zcoms));
        ddx = (pG->Phi[k][j][i+1] - pG->Phi[k][j][i-1]) / (2.0 * pG->dx1);
        ddy = (pG->Phi[k][j+1][i] - pG->Phi[k][j-1][i]) / (2.0 * pG->dx1);
        ddz = (pG->Phi[k+1][j][i] - pG->Phi[k-1][j][i]) / (2.0 * pG->dx1);
        locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = (ddx * (xd-xcoms) + ddy * (yd-ycoms) + ddz * (zd-zcoms)) / tv;
#ifdef RADIATION
    } else if (strcmp(expr,"Fedd")==0) {
        cc_pos(pG, i, j, k, &xd, &yd, &zd);
        ddx = (pG->Phi[k][j][i+1] - pG->Phi[k][j][i-1]) / (2.0 * pG->dx1);
        ddy = (pG->Phi[k][j+1][i] - pG->Phi[k][j-1][i]) / (2.0 * pG->dx1);
        ddz = (pG->Phi[k+1][j][i] - pG->Phi[k-1][j][i]) / (2.0 * pG->dx1);
        locarr[i-is+gis + gnx1*(j-js+gjs) + gnx1*gnx2*(k-ks+gks)] = (kappa_IR / c) * (pG->Urad[k][j][i].F1 * (xd-xcoms) + pG->Urad[k][j][i].F2 * (yd-ycoms) + pG->Urad[k][j][i].F3 * (zd-zcoms)) / (ddx * (xd-xcoms) + ddy * (yd-ycoms) + ddz * (zd-zcoms));
#endif
    } else {
	  ath_perr(-1,"Unknown data expression\n");
	}
      }
    }
  }   
        
#ifdef MPI_PARALLEL
  mpierr = MPI_Reduce(locarr, gridarr, gnx1*gnx2*gnx3, MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif

}

// Calculate a quantity evaluated in spherical shells around the centre of mass
static Real circpdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real rmeval, Real rmaxeval, int *gblsurfpdf, Real minsd, Real delsd, int npdf, int nang, char *expr, int massflag, int radflag)
{
    
    int i, j, k;
    long massupdate;
    int ntheta=nang,nphi=nang,nline=100;
    Real phimin, phimax, thetamin, thetamax, delphi, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp, delline;
    Real gradx, grady, gradz, rcoms, xcoms, ycoms, zcoms, tx, ty, tz, ndivs = 1.0e9;
    Real tempsd, tempden;
    Real totsum = 0.0, massline = 0.0, locmass;
    
    // Centre of mass of star particles
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Initialize the pdf to all zeroes
    for (i = 1; i <= npdf; i++) {
        gblsurfpdf[i-1] = 0;
    }
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            tempsd = 0.0;
            delline = rmaxeval / (nline - 1.0);
            massline = 0.0;
            for (k = 1; k <=nline; k++) {
                tx = gradx * (k - 1) * delline + xcoms;
                ty = grady * (k - 1) * delline + ycoms;
                tz = gradz * (k - 1) * delline + zcoms;
                tempden = intpgridpt(pG, gblarr, tx, ty, tz);
                locmass = tempden * delline * (k - 1) * delline * (k - 1) * delline * svalt * deltheta * delphi;
                if ((k - 1) * delline > rmeval) {
                    massline = massline + locmass;
                    if (radflag == 0) {
                        tempsd = tempsd + tempden * delline;
                    } else {
                        tempsd = tempsd + tempden * delline * locmass;
                    }
                }
            }
            if (radflag != 0) {
                tempsd = tempsd / massline;
            }
            
            if (strcmp(expr,"sd")==0) {
                totsum = totsum + tempsd * svalt * deltheta * delphi / (4.0 * 3.14159);
            } else if (strcmp(expr,"od")==0) {
                totsum = totsum + exp(-1.0 * tempsd * kappa_IR) * svalt * deltheta * delphi / (4.0 * 3.14159);
            } else {
                totsum = totsum + tempsd * svalt * deltheta * delphi / (4.0 * 3.14159);
            }
            
            // Update the outflow pdf according to this value
            if (massflag == 0) {
                intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfpdf, 1);
            } else {
                massupdate = (long)round(massline * ndivs / M_GMC);
                // ath_pout(0,"i = %d, j = %d, massupdate = %d, fupdate = %e, tempsd = %e\n", i, j, massupdate, massline * ndivs / M_GMC, tempsd);
                intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfpdf, massupdate);
            }
            
        }
    }
    
    return totsum;
}

// Calculate the amount of mass as a function of inclination angle
static Real anglemasspdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real rmeval, Real rmaxeval, int *gblanglepdf, int nang)
{
    
    int i, j, k;
    long massupdate;
    int ntheta=nang,nphi=nang,nline=100;
    Real phimin, phimax, thetamin, thetamax, delphi, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp, delline;
    Real gradx, grady, gradz, rcoms, xcoms, ycoms, zcoms, tx, ty, tz, ndivs = 1.0e9;
    Real tempden;
    Real totsum = 0.0, massline = 0.0, locmass;
    
    // Centre of mass of star particles
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Initialize the pdf to all zeroes
    for (i = 1; i <= ntheta; i++) {
        gblanglepdf[i-1] = 0;
    }
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            delline = rmaxeval / (nline - 1.0);
            massline = 0.0;
            for (k = 1; k <=nline; k++) {
                tx = gradx * (k - 1) * delline + xcoms;
                ty = grady * (k - 1) * delline + ycoms;
                tz = gradz * (k - 1) * delline + zcoms;
                tempden = intpgridpt(pG, gblarr, tx, ty, tz);
                locmass = tempden * delline * (k - 1) * delline * (k - 1) * delline * svalt * deltheta * delphi;
                if ((k - 1) * delline > rmeval) {
                    massline = massline + locmass;
                }
            }
            
            totsum = totsum + massline;
            
            // Update the outflow pdf according to this value
            massupdate = (long)round(massline * ndivs / M_GMC);
            // ath_pout(0,"i = %d, j = %d, massupdate = %d, fupdate = %e, tempsd = %e\n", i, j, massupdate, massline * ndivs / M_GMC, valspht);
            intupdarr(valspht, thetamin, delspht, ntheta, gblanglepdf, massupdate);
            
        }
    }
    // ath_error("Premature Exit\n");
    
    return totsum;
}

// Calculate a quantity evaluated in spherical shells around the centre of mass
static Real circtlimpdf(const GridS *pG, Real *gblarr, Real *xcomvec, Real *tlims, Real rmeval, Real rmaxeval, int *gblsurfpdf, Real minsd, Real delsd, int npdf, int nang, char *expr, int massflag, int radflag)
{
    
    int i, j, k;
    long massupdate;
    int ntheta=nang,nphi=nang,nline=100;
    Real phimin, phimax, thetamina, thetamaxa, thetaminb, thetamaxb, delphi, delsphta, delsphtb, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp, delline;
    Real gradx, grady, gradz, rcoms, xcoms, ycoms, zcoms, tx, ty, tz, ndivs = 1.0e9;
    Real tempsd, tempden;
    Real totsum = 0.0, massline = 0.0, locmass;
    
    // Centre of mass of star particles
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Initialize the pdf to all zeroes
    for (i = 1; i <= npdf; i++) {
        gblsurfpdf[i-1] = 0;
    }
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamina = 0.5 * (cos(tlims[0]) + 1);
    thetamaxa = 0.5 * (cos(tlims[1]) + 1);
    delsphta = 2.0 * (thetamaxa - thetamina) / (ntheta);
    thetaminb = 0.5 * (cos(tlims[2]) + 1);
    thetamaxb = 0.5 * (cos(tlims[3]) + 1);
    delsphtb = 2.0 * (thetamaxb - thetaminb) / (ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        if (i <= ntheta/2) {
            valspht = thetamina + (i - 1 + 0.5) * delsphta;
            delspht = delsphta;
        } else {
            valspht = thetaminb + (i - 1 + 0.5 - ntheta/2) * delsphtb;
            delspht = delsphtb;
        }
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            tempsd = 0.0;
            delline = rmaxeval / (nline - 1.0);
            massline = 0.0;
            for (k = 1; k <=nline; k++) {
                tx = gradx * (k - 1) * delline + xcoms;
                ty = grady * (k - 1) * delline + ycoms;
                tz = gradz * (k - 1) * delline + zcoms;
                tempden = intpgridpt(pG, gblarr, tx, ty, tz);
                locmass = tempden * delline * (k - 1) * delline * (k - 1) * delline * svalt * deltheta * delphi;
                if ((k - 1) * delline > rmeval) {
                    massline = massline + locmass;
                    if (radflag == 0) {
                        tempsd = tempsd + tempden * delline;
                    } else {
                        tempsd = tempsd + tempden * delline * locmass;
                    }
                }
            }
            if (radflag != 0) {
                tempsd = tempsd / massline;
            }
            
            if (strcmp(expr,"sd")==0) {
                totsum = totsum + tempsd * svalt * deltheta * delphi / (4.0 * 3.14159);
            } else if (strcmp(expr,"od")==0) {
                totsum = totsum + exp(-1.0 * tempsd * kappa_IR) * svalt * deltheta * delphi / (4.0 * 3.14159);
            } else {
                totsum = totsum + tempsd * svalt * deltheta * delphi / (4.0 * 3.14159);
            }
            
            // Update the outflow pdf according to this value
            if (massflag == 0) {
                intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfpdf, 1);
            } else {
                massupdate = (long)round(massline * ndivs / M_GMC);
                // ath_pout(0,"i = %d, j = %d, massupdate = %d, fupdate = %e, tempsd = %e\n", i, j, massupdate, massline * ndivs / M_GMC, tempsd);
                intupdarr(log10(tempsd), minsd, delsd, npdf, gblsurfpdf, massupdate);
            }
            
        }
    }
    
    return totsum;
}


// Calculate a quantity evaluated in spherical shells around the centre of mass
static void fltcircpdf(const GridS *pG, Real *gblarr, Real *augarr, Real *xcomvec, Real rmeval, Real rmaxeval, Real vlim, Real *gblvpdf, int *gblsurfpdf, Real minsd, Real delsd, Real minvel, Real delvel, int npdf, int nang, int massflag)
{

    int i, j, k;
    int ntheta=nang,nphi=nang,nline=100;
    long massupdate;
    Real phimin, phimax, thetamin, thetamax, delphi, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp, delline;
    Real gradx, grady, gradz, rcoms, xcoms, ycoms, zcoms, tx, ty, tz, ndivs = 1.0e9;
    Real tempsd, tempden, temptotden, temptotsd, tempv, temptotv;
    Real totsum = 0.0, massline = 0.0, locmass, locfac;

    // Centre of mass of star particles
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Initialize the pdf to all zeroes
    for (i = 1; i <= npdf; i++) {
        gblsurfpdf[i-1] = 0;
        gblvpdf[i-1] = 0.0;
    }
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            tempsd = 0.0;
            delline = rmaxeval / (nline - 1.0);
            temptotden = 0.0;
            temptotsd = 0.0;
            temptotv = 0.0;
            massline = 0.0;
            for (k = 1; k <=nline; k++) {
                tx = gradx * (k - 1) * delline + xcoms;
                ty = grady * (k - 1) * delline + ycoms;
                tz = gradz * (k - 1) * delline + zcoms;
                tempden = intpgridpt(pG, augarr, tx, ty, tz);
                tempsd = tempden * delline;
                locfac = delline * (k - 1) * delline * (k - 1) * delline * svalt * deltheta * delphi;
                locmass = tempden * locfac;
                tempv = intpgridpt(pG, gblarr, tx, ty, tz);
                if ((k - 1) * delline > rmeval) {
                    massline = massline + locmass;
                    temptotsd = temptotsd + tempsd;
                    temptotden = temptotden + locmass;
                    temptotv = temptotv + tempv * locfac;
                }
            }
            
            if (temptotden > 0.0) {
                temptotv = temptotv / temptotden;
            }
            
            // Update the outflow pdf according to this value
            if (temptotv > vlim) {
                if (massflag == 0) {
                    intupdarr(log10(temptotsd), minsd, delsd, npdf, gblsurfpdf, 1);
                    fltupdarr(log10(temptotsd), minsd, delsd, npdf, gblvpdf, temptotv);
                } else {
                    massupdate = (long)round(massline * ndivs / M_GMC);
                    intupdarr(log10(temptotsd), minsd, delsd, npdf, gblsurfpdf, massupdate);
                    fltupdarr(log10(temptotsd), minsd, delsd, npdf, gblvpdf, temptotv * massline * ndivs / M_GMC);
                }
            }
        }
    }

}

// Calculate the mean radius along sightlines in spherical shells around the centre of mass
static void fltradcircpdf(const GridS *pG, Real *augarr, Real *xcomvec, Real rmeval, Real rmaxeval, Real *gblvpdf, int *gblsurfpdf, Real minsd, Real delsd, Real minvel, Real delvel, int npdf, int nang, int massflag)
{
    
    int i, j, k;
    int ntheta=nang,nphi=nang,nline=100;
    long massupdate;
    Real phimin, phimax, thetamin, thetamax, delphi, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp, delline;
    Real gradx, grady, gradz, rcoms, xcoms, ycoms, zcoms, tx, ty, tz, ndivs = 1.0e9;
    Real tempsd, tempden, temptotden, temptotsd, tempv, temptotv;
    Real totsum = 0.0, massline = 0.0, locmass, locfac;
    
    // Centre of mass of star particles
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Initialize the pdf to all zeroes
    for (i = 1; i <= npdf; i++) {
        gblsurfpdf[i-1] = 0;
        gblvpdf[i-1] = 0.0;
    }
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            tempsd = 0.0;
            delline = rmaxeval / (nline - 1.0);
            temptotden = 0.0;
            temptotsd = 0.0;
            temptotv = 0.0;
            massline = 0.0;
            for (k = 1; k <=nline; k++) {
                tx = gradx * (k - 1) * delline + xcoms;
                ty = grady * (k - 1) * delline + ycoms;
                tz = gradz * (k - 1) * delline + zcoms;
                tempden = intpgridpt(pG, augarr, tx, ty, tz);
                tempsd = tempden * delline;
                locfac = delline * (k - 1) * delline * (k - 1) * delline * svalt * deltheta * delphi;
                locmass = tempden * locfac;
                if ((k - 1) * delline > rmeval) {
                    temptotsd = temptotsd + tempsd;
                    massline = massline + locmass;
                    temptotv = temptotv + locmass * (k - 1) * delline;
                }
            }
            
            if (massline > 0.0) {
                temptotv = temptotv / massline;
            }
            
            if (massflag == 0) {
                intupdarr(log10(temptotsd), minsd, delsd, npdf, gblsurfpdf, 1);
                fltupdarr(log10(temptotsd), minsd, delsd, npdf, gblvpdf, temptotv);
            } else {
                massupdate = (long)round(massline * ndivs / M_GMC);
                intupdarr(log10(temptotsd), minsd, delsd, npdf, gblsurfpdf, massupdate);
                fltupdarr(log10(temptotsd), minsd, delsd, npdf, gblvpdf, temptotv * massline * ndivs / M_GMC);
            }
            
        }
    }
    
}

// Calculate a field evaluated on a sphere for a set of given input radii
static void circfield(const GridS *pG, Real *gblarr, Real *rflux, Real *gblflux, Real *gblarea, Real *xcomvec, int nang, int nrs, int posflag)
{
    
    int i, j, k;
    int ntheta=nang,nphi=nang;
    Real phimin, phimax, thetamin, thetamax, delphi, delr, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp;
    Real gradx, grady, gradz, tx, ty, tz, txp, typ, tzp, txn, tyn, tzn, xcoms, ycoms, zcoms, rcoms;
    Real tempsd, tempden, rmeval, tempdenp, tempdenm;
    
    // Initialize the flux to all zeroes
    for (i = 1; i <= ntheta*nphi*nrs; i++) {
        gblflux[i-1] = 0.0;
        gblarea[i-1] = 0.0;
    }
    
    // Centre of mass
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    delr = pG->dx1;
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            
            for (k = 1; k <= nrs; k++) {
                rmeval = rflux[k - 1];
                tx = gradx * rmeval + xcoms;
                ty = grady * rmeval + ycoms;
                tz = gradz * rmeval + zcoms;
                tempden = intpgridpt(pG, gblarr, tx, ty, tz);
                gblarea[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = rmeval * rmeval * svalt * deltheta * delphi;
                
                if (posflag == 0) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden;
                } else if (posflag == 1) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * gradx * rmeval;
                } else if (posflag == 2) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * grady * rmeval;
                } else if (posflag == 3) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * gradz * rmeval;
                } else {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden;
                }
            }
            
        }
    }
}

// Calculate a field evaluated on a sphere for a set of given input radii on a grid corresponding to one processor
static void loclcircfield(const GridS *pG, Real *rflux, Real *gblflux, Real *gblarea, Real *xcomvec, int nang, int nrs, int posflag, char *expr)
{
    
    int i, j, k;
    int ntheta=nang,nphi=nang;
    Real phimin, phimax, thetamin, thetamax, delphi, delspht, valspht, valt, svalt, cvalt, deltheta, valp, svalp, cvalp;
    Real gradx, grady, gradz, tx, ty, tz, xcoms, ycoms, zcoms, rcoms;
    Real tempsd, tempden, rmeval;
    
    // Initialize the flux to all zeroes
    for (i = 1; i <= ntheta*nphi*nrs; i++) {
        gblflux[i-1] = 0.0;
        gblarea[i-1] = 0.0;
    }
    
    // Centre of mass
    xcoms = xcomvec[0];
    ycoms = xcomvec[1];
    zcoms = xcomvec[2];
    rcoms = sqrt(xcoms * xcoms + ycoms * ycoms + zcoms * zcoms);
    
    // Loop over successive values of theta and phi in bins of fixed solid angle
    phimin = 0.0;
    phimax = 2.0 * PI;
    thetamin = 1.0;
    thetamax = 0.0;
    delspht = (thetamax - thetamin) / (1.0 * ntheta);
    delphi = (phimax - phimin) / (1.0 * nphi);
    for (i = 1; i <= ntheta; i++) {
        valspht = thetamin + (i - 1 + 0.5) * delspht;
        valt = acos(2.0 * valspht - 1.0);
        svalt = sin(valt);
        cvalt = cos(valt);
        deltheta = -2.0 * delspht / sqrt(1.0 - (2.0 * valspht - 1.0) * (2.0 * valspht - 1.0));
        for (j = 1; j <= nphi; j++) {
            valp = phimin + (j - 1 + 0.5) * delphi;
            svalp = sin(valp);
            cvalp = cos(valp);
            gradx = svalt * cvalp;
            grady = svalt * svalp;
            gradz = cvalt;
            
            for (k = 1; k <= nrs; k++) {
                rmeval = rflux[k - 1];
                tx = gradx * rmeval + xcoms;
                ty = grady * rmeval + ycoms;
                tz = gradz * rmeval + zcoms;
                tempden = intploclgridpt(pG, tx, ty, tz, expr);
                gblarea[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = rmeval * rmeval * svalt * deltheta * delphi;
                
                if (posflag == 0) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden;
                } else if (posflag == 1) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * gradx * rmeval;
                } else if (posflag == 2) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * grady * rmeval;
                } else if (posflag == 3) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * gradz * rmeval;
                } else if (posflag == 4) {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden * gradz * rmeval;
                } else {
                    gblflux[(k - 1) * nphi * ntheta + (i - 1) * nphi + j - 1] = tempden;
                }
            }
            
        }
    }
}

// Interpolate a grid point at location (tx, ty, tz)
static Real intpgridpt(const GridS *pG, Real *gridarr, Real tx, Real ty, Real tz)
{

  int ilwr, jlwr, klwr;
  Real txl, tyl, tzl, xd, yd, zd, trip[7];
    
  ilwr = (int)(pG->is + floor((tx - pG->MinX[0])/pG->dx1 - 0.5));
  jlwr = (int)(pG->js + floor((ty - pG->MinX[1])/pG->dx2 - 0.5));
  klwr = (int)(pG->ks + floor((tz - pG->MinX[2])/pG->dx3 - 0.5));
  txl = pG->MinX[0] + ((Real)(ilwr - pG->is) + 0.5)*pG->dx1;
  tyl = pG->MinX[1] + ((Real)(jlwr - pG->js) + 0.5)*pG->dx2;
  tzl = pG->MinX[2] + ((Real)(klwr - pG->ks) + 0.5)*pG->dx3;
  xd = (tx - txl) / pG->dx1;
  yd = (ty - tyl) / pG->dx2;
  zd = (tz - tzl) / pG->dx3;

  // Calculate the density at this interpolated point
  trip[0] = gridarr[ilwr-is+gis + gnx1*(jlwr-js+gjs) + gnx1*gnx2*(klwr-ks+gks)] * (1.0 - xd) + gridarr[ilwr-is+gis+1 + gnx1*(jlwr-js+gjs) + gnx1*gnx2*(klwr-ks+gks)] * xd;
  trip[1] = gridarr[ilwr-is+gis + gnx1*(jlwr-js+gjs+1) + gnx1*gnx2*(klwr-ks+gks)] * (1.0 - xd) + gridarr[ilwr-is+gis+1 + gnx1*(jlwr-js+gjs+1) + gnx1*gnx2*(klwr-ks+gks)] * xd;
  trip[2] = gridarr[ilwr-is+gis + gnx1*(jlwr-js+gjs) + gnx1*gnx2*(klwr-ks+gks+1)] * (1.0 - xd) + gridarr[ilwr-is+gis+1 + gnx1*(jlwr-js+gjs) + gnx1*gnx2*(klwr-ks+gks+1)] * xd;
  trip[3] = gridarr[ilwr-is+gis + gnx1*(jlwr-js+gjs+1) + gnx1*gnx2*(klwr-ks+gks+1)] * (1.0 - xd) + gridarr[ilwr-is+gis+1 + gnx1*(jlwr-js+gjs+1) + gnx1*gnx2*(klwr-ks+gks+1)] * xd;
  trip[4] = trip[0] * (1.0 - yd) + trip[1] * yd;
  trip[5] = trip[2] * (1.0 - yd) + trip[3] * yd;
  trip[6] = trip[4] * (1.0 - zd) + trip[5] * zd;

  return trip[6];
}

// Defunct function
static Real intploclgridpt(const GridS *pG, Real tx, Real ty, Real tz, char *expr)
{
    
    int ilwr, jlwr, klwr;
    Real txl, tyl, tzl, xd, yd, zd, trip[7];
    
    if (tx < pG->MinX[0] || tx > pG->MaxX[0] || ty < pG->MinX[1] || ty > pG->MaxX[1] || tz < pG->MinX[2] || tz > pG->MaxX[2]) {
        return 0.0;
    } else {
        ilwr = (int)(pG->is + floor((tx - pG->MinX[0])/pG->dx1 - 0.5));
        jlwr = (int)(pG->js + floor((ty - pG->MinX[1])/pG->dx2 - 0.5));
        klwr = (int)(pG->ks + floor((tz - pG->MinX[2])/pG->dx3 - 0.5));
        txl = pG->MinX[0] + ((Real)(ilwr - pG->is) + 0.5)*pG->dx1;
        tyl = pG->MinX[1] + ((Real)(jlwr - pG->js) + 0.5)*pG->dx2;
        tzl = pG->MinX[2] + ((Real)(klwr - pG->ks) + 0.5)*pG->dx3;
        xd = (tx - txl) / pG->dx1;
        yd = (ty - tyl) / pG->dx2;
        zd = (tz - tzl) / pG->dx3;
        
        // Calculate the density at this interpolated point
        if (strcmp(expr,"d")==0) {
            trip[0] = pG->U[klwr][jlwr][ilwr].d * (1.0 - xd) + pG->U[klwr][jlwr][ilwr+1].d * xd;
            trip[1] = pG->U[klwr][jlwr+1][ilwr].d * (1.0 - xd) + pG->U[klwr][jlwr+1][ilwr+1].d * xd;
            trip[2] = pG->U[klwr+1][jlwr][ilwr].d * (1.0 - xd) + pG->U[klwr+1][jlwr][ilwr+1].d * xd;
            trip[3] = pG->U[klwr+1][jlwr+1][ilwr].d * (1.0 - xd) + pG->U[klwr+1][jlwr+1][ilwr+1].d * xd;
        } else if (strcmp(expr,"M1")==0) {
            trip[0] = pG->U[klwr][jlwr][ilwr].M1 * (1.0 - xd) + pG->U[klwr][jlwr][ilwr+1].M1 * xd;
            trip[1] = pG->U[klwr][jlwr+1][ilwr].M1 * (1.0 - xd) + pG->U[klwr][jlwr+1][ilwr+1].M1 * xd;
            trip[2] = pG->U[klwr+1][jlwr][ilwr].M1 * (1.0 - xd) + pG->U[klwr+1][jlwr][ilwr+1].M1 * xd;
            trip[3] = pG->U[klwr+1][jlwr+1][ilwr].M1 * (1.0 - xd) + pG->U[klwr+1][jlwr+1][ilwr+1].M1 * xd;
        } else if (strcmp(expr,"M2")==0) {
            trip[0] = pG->U[klwr][jlwr][ilwr].M2 * (1.0 - xd) + pG->U[klwr][jlwr][ilwr+1].M2 * xd;
            trip[1] = pG->U[klwr][jlwr+1][ilwr].M2 * (1.0 - xd) + pG->U[klwr][jlwr+1][ilwr+1].M2 * xd;
            trip[2] = pG->U[klwr+1][jlwr][ilwr].M2 * (1.0 - xd) + pG->U[klwr+1][jlwr][ilwr+1].M2 * xd;
            trip[3] = pG->U[klwr+1][jlwr+1][ilwr].M2 * (1.0 - xd) + pG->U[klwr+1][jlwr+1][ilwr+1].M2 * xd;
        } else if (strcmp(expr,"M3")==0) {
            trip[0] = pG->U[klwr][jlwr][ilwr].M3 * (1.0 - xd) + pG->U[klwr][jlwr][ilwr+1].M3 * xd;
            trip[1] = pG->U[klwr][jlwr+1][ilwr].M3 * (1.0 - xd) + pG->U[klwr][jlwr+1][ilwr+1].M3 * xd;
            trip[2] = pG->U[klwr+1][jlwr][ilwr].M3 * (1.0 - xd) + pG->U[klwr+1][jlwr][ilwr+1].M3 * xd;
            trip[3] = pG->U[klwr+1][jlwr+1][ilwr].M3 * (1.0 - xd) + pG->U[klwr+1][jlwr+1][ilwr+1].M3 * xd;
#ifdef RADIATION
        } else if (strcmp(expr,"F1")==0) {
            trip[0] = pG->Urad[klwr][jlwr][ilwr].F1 * (1.0 - xd) + pG->Urad[klwr][jlwr][ilwr+1].F1 * xd;
            trip[1] = pG->Urad[klwr][jlwr+1][ilwr].F1 * (1.0 - xd) + pG->Urad[klwr][jlwr+1][ilwr+1].F1 * xd;
            trip[2] = pG->Urad[klwr+1][jlwr][ilwr].F1 * (1.0 - xd) + pG->Urad[klwr+1][jlwr][ilwr+1].F1 * xd;
            trip[3] = pG->Urad[klwr+1][jlwr+1][ilwr].F1 * (1.0 - xd) + pG->Urad[klwr+1][jlwr+1][ilwr+1].F1 * xd;
        } else if (strcmp(expr,"F2")==0) {
            trip[0] = pG->Urad[klwr][jlwr][ilwr].F2 * (1.0 - xd) + pG->Urad[klwr][jlwr][ilwr+1].F2 * xd;
            trip[1] = pG->Urad[klwr][jlwr+1][ilwr].F2 * (1.0 - xd) + pG->Urad[klwr][jlwr+1][ilwr+1].F2 * xd;
            trip[2] = pG->Urad[klwr+1][jlwr][ilwr].F2 * (1.0 - xd) + pG->Urad[klwr+1][jlwr][ilwr+1].F2 * xd;
            trip[3] = pG->Urad[klwr+1][jlwr+1][ilwr].F2 * (1.0 - xd) + pG->Urad[klwr+1][jlwr+1][ilwr+1].F2 * xd;
        } else if (strcmp(expr,"F3")==0) {
            trip[0] = pG->Urad[klwr][jlwr][ilwr].F3 * (1.0 - xd) + pG->Urad[klwr][jlwr][ilwr+1].F3 * xd;
            trip[1] = pG->Urad[klwr][jlwr+1][ilwr].F3 * (1.0 - xd) + pG->Urad[klwr][jlwr+1][ilwr+1].F3 * xd;
            trip[2] = pG->Urad[klwr+1][jlwr][ilwr].F3 * (1.0 - xd) + pG->Urad[klwr+1][jlwr][ilwr+1].F3 * xd;
            trip[3] = pG->Urad[klwr+1][jlwr+1][ilwr].F3 * (1.0 - xd) + pG->Urad[klwr+1][jlwr+1][ilwr+1].F3 * xd;
#endif
        } else if (strcmp(expr,"Phi")==0) {
            trip[0] = pG->Phi[klwr][jlwr][ilwr] * (1.0 - xd) + pG->Phi[klwr][jlwr][ilwr+1] * xd;
            trip[1] = pG->Phi[klwr][jlwr+1][ilwr] * (1.0 - xd) + pG->Phi[klwr][jlwr+1][ilwr+1] * xd;
            trip[2] = pG->Phi[klwr+1][jlwr][ilwr] * (1.0 - xd) + pG->Phi[klwr+1][jlwr][ilwr+1] * xd;
            trip[3] = pG->Phi[klwr+1][jlwr+1][ilwr] * (1.0 - xd) + pG->Phi[klwr+1][jlwr+1][ilwr+1] * xd;
        } else {
            return 0.0;
        }
        trip[4] = trip[0] * (1.0 - yd) + trip[1] * yd;
        trip[5] = trip[2] * (1.0 - yd) + trip[3] * yd;
        trip[6] = trip[4] * (1.0 - zd) + trip[5] * zd;
        
        return trip[6];
    }
}

// Given an input density, min-density, max-density, #of density bins in PDF, output array for PDF and integer value with which
// to update PDF at this density, augments datarr by this value
static void intupdarr(Real tval, Real minval, Real delval, int nvals, int *datarr, int updval)
{

  int iden;

  iden = (int) ((tval - minval) / delval);
  if (iden >= nvals) {
    iden=nvals-1;
  }
  if (iden < 0) {
    iden=0;
  }
  // ath_pout("tval = (%e, %e, %e), iden = %d\n", minval, tval, delval, iden);
  datarr[iden]=datarr[iden]+updval;

}

// Same but for floating point arrays
static void fltupdarr(Real tval, Real minval, Real delval, int nvals, Real *datarr, Real updval)
{

  int iden;

  iden = (int) ((tval - minval) / delval);
  if (iden >= nvals) {
    iden=nvals-1;
  }
  if (iden < 0) {
    iden=0;
  }
  datarr[iden]=datarr[iden]+updval;

}


/*------------------------------------------------------------------------------
 *  TURBULENCE FUNCTIONS
 *----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
 *  Function pspect -- computes component of velocity with specific
 *  power spectrum in Fourier space determined by ispect
 *
 *  Velocity power spectrum returned in ampl
 *    klow   = multiple of 2 pi/L for cut-off at low  wavenumbers
 *    khigh  = multiple of 2 pi/L for cut-off at high wavenumbers
 *    expo   = exponent of power law
 *    ispect = integer flag which specifies spectrum
 *  Note that the fourier amplitudes are stored in an array with no
 *  ghost zones
 */
static void pspect(ath_fft_data *ampl)
{
  int i,j,k;
  Real q1,q2,q3,q4,q5;
  
  for (k=0; k<gnx3; k++) {
    for (j=0; j<gnx2; j++) {
      for (i=0; i<gnx1; i++) {
        /* Calculate the wave vector magnitude, k, in global coordinates,
         * even if (i-gis,j-gjs,k-gks) is not on the local Grid. */
        q4 = KWVM(i-gis,j-gjs,k-gks);
        /* If k is within the cutoff range, generate 3 random numbers. 
         * This is done so that different processor numbers/topologies will
         * produce the exact same random number sequences. */
        if ((q4 > klow) && (q4 < khigh)) {
          q1 = ran2(&rseed);
          q2 = ran2(&rseed);
          q3 = ran2(&rseed);
        }
        
        /* If (i-gis,j-gjs,k-gks) is on the local Grid (in global coordinates), 
         * assign a wave amplitude (possibly 0). */
        if ((i>=gis) && (i<=gie) &&
            (j>=gjs) && (j<=gje) &&
            (k>=gks) && (k<=gke)) {
          /* If k is within the cutoff range, compute an amplitude. */
          if ((q4 > klow) && (q4 < khigh)) {
            q5 = sqrt(-2.0*log(q1+1.0e-20))*cos(2.0*PI*q2);
            ampl[OFST(i-gis,j-gjs,k-gks)][0] = q5*cos(2.0*PI*q3);
            ampl[OFST(i-gis,j-gjs,k-gks)][1] = q5*sin(2.0*PI*q3);
          
            /* Set power spectrum
             *   ispect=1: power law - original form
             *   ispect=2: form from Gammie&Ostriker
             */
            q4 *= dkx; /* Multiply by 2 pi/L */
            if (ispect == 1) {
              /* Decreasing power law */
              ampl[OFST(i-gis,j-gjs,k-gks)][0] /= pow(q4,(expo+2.0)/2.0);
              ampl[OFST(i-gis,j-gjs,k-gks)][1] /= pow(q4,(expo+2.0)/2.0);
            } else if (ispect == 2) {
              /* G&O form */
              ampl[OFST(i-gis,j-gjs,k-gks)][0] *= pow(q4,3.0)*exp(-4.0*q4/kpeak);
              ampl[OFST(i-gis,j-gjs,k-gks)][1] *= pow(q4,3.0)*exp(-4.0*q4/kpeak);
            }
          } else {
            /* Otherwise, introduce cut-offs at klow and khigh */
            ampl[OFST(i-gis,j-gjs,k-gks)][0] = 0.0;
            ampl[OFST(i-gis,j-gjs,k-gks)][1] = 0.0;
          }
        }
      }
    }
  }
  if (gis==0 && gie==0 && gjs==0 && gje==0 && gks==0 && gke==0) {
    ampl[0][0] = 0.0;
    ampl[0][1] = 0.0;
  }
  
  return;
}

/*------------------------------------------------------------------------------
 *  Function project
 *
 *  Makes velocity perturbations divergence free
 */
static void project()
{
  int i,j,k,m,ind;
  double kap[3], kapn[3], mag;
  ath_fft_data dot;
  
  /* Project off non-solenoidal component of velocity */
  for (k=0; k<nx3; k++) {
    kap[2] = sin(2.0*PI*(gks+k)/gnx3);
    for (j=0; j<nx2; j++) {
      kap[1] = sin(2.0*PI*(gjs+j)/gnx2);
      for (i=0; i<nx1; i++) {
        if (((gis+i)+(gjs+j)+(gks+k)) != 0) {
          kap[0] = sin(2.0*PI*(gis+i)/gnx1);
          ind = OFST(i,j,k);
          
          /* make kapn a unit vector */
          mag = sqrt(SQR(kap[0]) + SQR(kap[1]) + SQR(kap[2]));
          for (m=0; m<3; m++) kapn[m] = kap[m] / mag;
          
          /* find fv_0 dot kapn */
          dot[0] = fv1[ind][0]*kapn[0]+fv2[ind][0]*kapn[1]+fv3[ind][0]*kapn[2];
          dot[1] = fv1[ind][1]*kapn[0]+fv2[ind][1]*kapn[1]+fv3[ind][1]*kapn[2];
          
          /* fv = fv_0 - (fv_0 dot kapn) * kapn */
          fv1[ind][0] -= dot[0]*kapn[0];
          fv2[ind][0] -= dot[0]*kapn[1];
          fv3[ind][0] -= dot[0]*kapn[2];
          
          fv1[ind][1] -= dot[1]*kapn[0];
          fv2[ind][1] -= dot[1]*kapn[1];
          fv3[ind][1] -= dot[1]*kapn[2];
        }
      }
    }
  }
  
  return;
}

/*------------------------------------------------------------------------------
 *  Function transform
 *
 *  Generate velocities from fourier transform
 */
static inline void transform()
{
  /* Transform velocities from k space to physical space */
  ath_3d_fft(plan, fv1);
  ath_3d_fft(plan, fv2);
  ath_3d_fft(plan, fv3);
  
  /* Should technically renormalize (divide by gnx1*gnx2*gnx3) here, but
   * since we're going to renormalize to get the desired energy injection
   * rate anyway, there's no point */
  
  return;
}

/*------------------------------------------------------------------------------
 *  Function generate
 *
 *  Generate the velocity perturbations
 */
static inline void generate()
{
  /* Generate new perturbations following appropriate power spectrum */
  pspect(fv1);
  pspect(fv2);
  pspect(fv3);
  
  /* Require div V = 0 */
  //  project();
  
  /* Transform perturbations to real space, but don't normalize until
   * just before we apply them in perturb() */
  transform();
  
  return;
}

/*------------------------------------------------------------------------------
 *  Function perturb
 *
 *  Shifts velocities so no net momentum change, normalizes to keep
 *  dedt fixed, and then sets velocities (momenta)
 */
static void perturb(GridS *pG, Real dt)
{
  int i,j,k;
  int ind, mpierr;
  Real dvol, s, de, qa, v1, v2, v3;
  Real t0,t1,t2,t3,aa,bb,cc;
  Real m[4], gm[4];
  
  /* Set the velocities in real space */
  dvol = 1.0/((Real)(gnx1*gnx2*gnx3));
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        ind = OFST(i-is,j-js,k-ks);
        dv1[k][j][i] = fv1[ind][0]*dvol;
        dv2[k][j][i] = fv2[ind][0]*dvol;
        dv3[k][j][i] = fv3[ind][0]*dvol;
      }
    }
  }
  
  /* Calculate net momentum pertubation components t1, t2, t3 */
  t0 = 0.0;  t1 = 0.0;  t2 = 0.0;  t3 = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        t0 += pG->U[k][j][i].d;
        
        /* The net momentum perturbation */
        t1 += pG->U[k][j][i].d * dv1[k][j][i];
        t2 += pG->U[k][j][i].d * dv2[k][j][i];
        t3 += pG->U[k][j][i].d * dv3[k][j][i];
      }
    }
  }
  
#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = t0;  m[1] = t1;  m[2] = t2;  m[3] = t3;
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[radpargrav]: MPI_Allreduce error = %d\n", mpierr);
  t0 = gm[0];  t1 = gm[1];  t2 = gm[2];  t3 = gm[3];
#endif /* MPI_PARALLEL */
  
  /* Subtract the mean velocity perturbation so that the net momentum
   * perturbation is zero. */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        dv1[k][j][i] -= t1/t0;
        dv2[k][j][i] -= t2/t0;
        dv3[k][j][i] -= t3/t0;
      }
    }
  }
  
  /* Calculate unscaled energy of perturbations */
  aa = 0.0;  bb = 0.0;  cc = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        /* Calculate velocity pertubation at cell center from
         * perturbations at cell faces */
        v1 = dv1[k][j][i];
        v2 = dv2[k][j][i];
        v3 = dv3[k][j][i];
        
        aa += 0.5*(pG->U[k][j][i].d)*(SQR(v1) + SQR(v2) + SQR(v3));
        bb += (pG->U[k][j][i].M1)*v1 + (pG->U[k][j][i].M2)*v2 +
	  (pG->U[k][j][i].M3)*v3;
        cc += 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) +
                   SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
      }
    }
  }
  
#ifdef MPI_PARALLEL
  /* Sum the perturbations over all processors */
  m[0] = aa;  m[1] = bb;  m[2] = cc;
  mpierr = MPI_Allreduce(m, gm, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) ath_error("[radpargrav]: MPI_Allreduce error = %d\n", mpierr);
  aa = gm[0];  bb = gm[1];  cc = gm[2];
#endif /* MPI_PARALLEL */
  
  /* Rescale to give the correct energy injection rate */
  dvol = pG->dx1*pG->dx2*pG->dx3;
  aa = MAX(aa,1.0e-20);
  if (idrive == 0) {
    /* driven turbulence */
    cc = -dedt*dt/dvol;
  } else {
    /* decaying turbulence (all in one shot) */
    /* NOTE:  In this case, dedt is really the desired total kinetic energy */
    cc -= dedt/dvol;
  }
  s = (-bb + sqrt(SQR(bb) - 4.0*aa*cc))/(2.0*aa);
  if (isnan(s)) ath_error("[radpargrav]: s is NaN!\n");
  
  /* Apply momentum pertubations */
  t1 = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        qa = s*pG->U[k][j][i].d;
        pG->U[k][j][i].M1 += qa*dv1[k][j][i];
        pG->U[k][j][i].M2 += qa*dv2[k][j][i];
        pG->U[k][j][i].M3 += qa*dv3[k][j][i];
      }
    }
  }

  return;
}

/*------------------------------------------------------------------------------
 *  Function initialize
 *
 *  Allocate memory and initialize FFT plans
 *----------------------------------------------------------------------------*/
static void initialize(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int i,j,k;
  int nbuf, mpierr;
  Real kwv, kpara, kperp;
  char donedrive = 0;
  Real pc_cgs=3.085678e18,kms_cgs=1.0e5,mH_cgs=1.6733e-24;
  
  /*----------------------------------------------------------------------------
   * Variables within this block are stored globally, and used
   * within preprocessor macros.  Don't create variables with
   * these names within your function if you are going to use
   * OFST(), KCOMP(), or KWVM() within the function! */
  
  /* Get local grid size */
  nx1 = pG->Nx[0];
  nx2 = pG->Nx[1];
  nx3 = pG->Nx[2];
  
  /* Get extents of local grid */
  is = pG->is;  ie = pG->ie;  
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;
  
  /* Get global grid size */
  gnx1 = pD->Nx[0];
  gnx2 = pD->Nx[1];
  gnx3 = pD->Nx[2];
  
  /* Get extents of local FFT grid in global coordinates */
  gis = pG->Disp[0];  gie = pG->Disp[0]+pG->Nx[0]-1;
  gjs = pG->Disp[1];  gje = pG->Disp[1]+pG->Nx[1]-1;
  gks = pG->Disp[2];  gke = pG->Disp[2]+pG->Nx[2]-1;
  
  /* Get size of arrays with ghost cells */
  il = is-nghost*(pG->Nx[0]>1);  iu = ie+nghost*(pG->Nx[0]>1);  nx1gh = iu-il+1;
  jl = js-nghost*(pG->Nx[1]>1);  ju = je+nghost*(pG->Nx[1]>1);  nx2gh = ju-jl+1;
  kl = ks-nghost*(pG->Nx[2]>1);  ku = ke+nghost*(pG->Nx[2]>1);  nx3gh = ku-kl+1;

  /* Get spatial extents of local grid */
  cc_pos(pG,il,jl,kl,&x1min,&x2min,&x3min);
  cc_pos(pG,iu,ju,ku,&x1max,&x2max,&x3max);
 
  /* Get local grid length, volume */
  dx = MAX(MAX(pG->dx1,pG->dx2),pG->dx3);
  dV = pG->dx1*pG->dx2*pG->dx3;
  
  //  printf("proc %d:  nx1=%d, nx2=%d, nx3=%d\n", myID_Comm_world,nx1,nx2,nx3);
  //  printf("proc %d:  gis=%d, gjs=%d, gks=%d\n", myID_Comm_world,gis,gjs,gks);
  //  ath_pout(0,"gnx1=%d, gnx2=%d, gnx3=%d\n", gnx1,gnx2,gnx3);
  //  ath_pout(0,"nx1gh=%d, nx2gh=%d, nx3gh=%d\n", nx1gh,nx2gh,nx3gh);

  /*--------------------------------------------------------------------------*/
  /* Get input parameters */
  
  /* Parse input file */
  rho_small   = par_getd("problem", "rho_small");
  rho_ffac    = par_getd("problem", "rho_ffac");
  rho_fdt     = par_getd("problem", "rho_fdt");
  fluxrad_out = par_getd("problem", "fluxrad_out");
  surfd_out   = par_getd("problem", "surfd_out");
  M_GMC       = par_getd("problem", "M_GMC");
  rcloud      = par_getd("problem", "rcloud");
#ifdef RADIATION
  kappa_IR    = par_getd("problem", "kappa_IR");
  //  nradsrc     = par_geti("problem", "nradsrc");
  Psi         = par_getd("problem", "Psi");
#endif
  
  /* Initialize units
   *
   * This assumes inputs are converted such that code units are:
   * length unit = pc,
   * velocity unit = km/s
   * density unit = 1.4 m_H * 1 /cm^3  [ i.e. # density of H is in cm^-3]
   *
   * Thus:
   *
   * code length unit = 3.0856 e18 cm = 1 pc
   * code time unit   = pc/(km s^-1)=3.0856 e13 s = 0.978 Myr
   * code mass unit   = 1.4 mH * (pc/cm)^3 = 6.88 e31 g = 0.035 Msun
   */
  UnitS units;
  
  /* Code length unit */
  units.Lcode = pc_cgs;
  ath_pout(0,"length unit [cm]:  %e\n",units.Lcode);
  
  /* Code mass unit */
  units.Mcode = 1.4*mH_cgs*CUBE(pc_cgs);
  ath_pout(0,"mass unit [g]:  %e\n",units.Mcode);
  
  /* Code velocity unit */
  units.Vcode = kms_cgs;
  ath_pout(0,"velocity unit [cm s^-1]:  %e\n",units.Vcode);
  
  /* Code time unit */
  units.Tcode = units.Lcode/units.Vcode;
  ath_pout(0,"time unit [s]:  %e\n",units.Tcode);
  
  init_units(&units);
  ath_pout(1,"units.cm = %e\n",units.cm);
  ath_pout(1,"units.g  = %e\n",units.g);
  ath_pout(1,"units.s  = %e\n",units.s);

  Lx = (pD->MaxX[0] - pD->MinX[0])*pc_cgs*units.cm;
  Ly = (pD->MaxX[1] - pD->MinX[1])*pc_cgs*units.cm;
  Lz = (pD->MaxX[2] - pD->MinX[2])*pc_cgs*units.cm;

  /* Code density unit */
  ath_pout(0,"density unit [g cm^-3]:  %e\n",units.Dcode);
  
  /* Code energy unit */
  ath_pout(0,"1 erg in code units = %e\n",units.erg);
  
  /* Code temperature unit */
  /* Chosen so that aR=1 in code units */
  ath_pout(0,"1 K in code units = %e\n",units.K);

#ifdef RADIATION
  /* Scale kappa to code units */
  kappa_IR *= SQR(units.cm)/units.g;
  ath_pout(0,"kappa_IR in code units = %e\n",kappa_IR);
  ath_pout(0,"source emission is resolved over %e zones\n",rsrc/dx);
#endif
  
  /* Scale M_GMC to code units */
  M_GMC *= units.Msun;
  ath_pout(0,"M_GMC in code units = %e\n",M_GMC);
  
  /* Scale rcloud to code units */
  rcloud *= units.pc;
  ath_pout(0,"rcloud in code units = %e\n",rcloud);  
  
  rho_cloud = rho(rcloud);
  ath_pout(0,"rho_cloud in code units = %e\n",rho_cloud);
  rho_small *= rho(rcloud);  /* scale rho_small by density maximum */
  ath_pout(0,"rho_small in code units = %e\n",rho_small);
  
#ifdef RADIATION
  /* Set gas constant in code units */
  Rgas = units.kB/(1.4*units.mH);
  ath_pout(0,"Rgas in code units = %e\n",Rgas);
  
  /* Set radiation constant in code units */
  /* aR = aR_cgs [erg cm^-3 K^-4] */
  aR = units.aR;
  ath_pout(0,"aR in code units = %e\n",aR);
  
  /* Set speed of light in code units */
  c = units.c;
  ath_pout(0,"c in code units = %e\n",c);
  
  /* Compute optical depth across cloud */
  tau_cloud = rho_cloud*kappa_IR*2.0*rcloud;
  ath_pout(0,"cloud optical depth ~ %e\n",tau_cloud);
#endif
  
  /* Scale Psi to code units */
  Psi *= units.erg/(units.s*units.g);
  ath_pout(0,"Psi in code units = %e\n",Psi);
    
#ifdef SELF_GRAVITY
  /* Set gravity constant */
  /* G = G_cgs [cm^3 g^-1 s^-2] */
  four_pi_G = 4.0*PI*units.G;
  grav_mean_rho = 0.0;
  ath_pout(0,"4*pi*G in code units = %e\n",four_pi_G);
  
  ath_pout(0,"rho_Truelove in code units = %e\n",
           SQR(PI)*Iso_csound2/(4.0*four_pi_G*SQR(dx)));
  ath_pout(0,"rho_LP in code units = %e\n",
           8.86*Iso_csound2*4.0/(four_pi_G*SQR(dx)));
  /* Ensure the Jeans length is resolved by at least 4 zones */
  L_Jeans = Iso_csound*2.0*PI/sqrt(four_pi_G*rho_cloud);
  ath_pout(0,"L_Jeans in code units ~ %e, resolved over %1.1f zones\n",
           L_Jeans,L_Jeans/dx);
  if (L_Jeans/dx < 0.04)
    ath_error("[radpargrav]:  The Jeans length is not resolved!\n");
  ath_pout(0,"t_ff in code units = %e\n",sqrt(3.0*PI*PI/(8.0*four_pi_G*rho_cloud)));
#endif
  
  
  /*--------------------------------------------------------------------------*/
  /* Set up perturbation parameters */
  
  /* interval for generating new driving spectrum; also interval for
   * driving when IMPULSIVE_DRIVING is used */
  dtdrive = par_getd("problem","dtdrive");
#ifdef MHD
  /* magnetic field strength */
  beta = par_getd("problem","beta");
  /* beta = isothermal pressure/magnetic pressure */
  // B0 = sqrt(2*Iso_csound2*rho_cloud/beta);
  B0 = sqrt(8.0*PI*Iso_csound2*rho_cloud/beta);
#endif /* MHD */
  /* energy injection rate */
  //  dedt = par_getd("problem","dedt");
  //  dedt *= units.erg;
  /* Set dedt(=Ekin,tot) such that Ekin = 0.5*alpha_vir*Egrav, where
   * Egrav = (3/5)*G*M^2/R  */
  Egrav = 0.6*(four_pi_G/(4.0*PI))*SQR(M_GMC)/rcloud;
  dedt = 0.5*par_getd("problem","alpha_vir")*Egrav;
  ath_pout(0,"Ekin in code units = %e\n",dedt);
  ath_pout(0,"Egrav in code units = %e\n",Egrav);
  v_turb = sqrt(2.0*dedt/M_GMC);
  ath_pout(0,"v_turb in code units = %e\n",v_turb);
#ifdef RADIATION
  ath_pout(0,"crad >> %e km s^-1 is recommended!\n",tau_cloud*(v_turb+Iso_csound));
#endif
  
  /* parameters for spectrum */
  ispect = par_geti("problem","ispect");
  if (ispect == 1) {
    expo = par_getd("problem","expo");
  } else if (ispect == 2) {
    kpeak = par_getd("problem","kpeak")*2.0*PI;
  } else {
    ath_error("[radpargrav]:  Invalid value for ispect\n");
  }
  
  /* Cutoff wavenumbers of spectrum */
  klow = par_getd("problem","klow"); /* in integer units */
  khigh = par_getd("problem","khigh"); /* in integer units */
  //  dkx = 2.0*PI/(pG->dx1*gnx1); /* convert k from integer */
  dkx = 2.0*PI/MAX(MAX(pG->dx1*gnx1,pG->dx2*gnx2),pG->dx3*gnx3); /* convert k from integer */
  
  /* Driven or decaying */
  idrive = par_geti("problem","idrive");
  if ((idrive < 0) || (idrive > 1))
    ath_error("[radpargrav]:  Invalid value for idrive\n");
  /* If restarting with decaying turbulence, no driving necessary. */
  if ((idrive == 1) && (pG->time > 0)) {
    donedrive = 1;
  }
  
  ath_pout(0,"Allocating memory for velocities\n");
  if (donedrive == 0) {
    /* Allocate memory for components of velocity perturbation */
    if ((dv1=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[radpargrav]: Error allocating memory for vel pert\n");
    }
    if ((dv2=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[radpargrav]: Error allocating memory for vel pert\n");
    }
    if ((dv3=(Real***)calloc_3d_array(nx3gh,nx2gh,nx1gh,sizeof(Real)))==NULL) {
      ath_error("[radpargrav]: Error allocating memory for vel pert\n");
    }
  }
  
  /* Initialize the FFT plan */
  ath_pout(0,"Initializing FFT Plan\n");
  plan = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
  
  /* Allocate memory for FFTs */
  ath_pout(0,"Allocating memory for FFTs\n");
  if (donedrive == 0) {
    fv1 = ath_3d_fft_malloc(plan);
    fv2 = ath_3d_fft_malloc(plan);
    fv3 = ath_3d_fft_malloc(plan);
  }
  ath_pout(0,"Finished Allocating memory for FFTs\n");
  
  /* Enroll outputs */
  dump_history_enroll(hst_mass_stars,"<mass_stars>");
#ifdef RADIATION
  dump_history_enroll(hst_F_out,"<F_out>");
#endif
  dump_history_enroll(hst_Egrav_gas,"<dEgrav_gas>");
  dump_history_enroll(hst_Egrav_stars,"<dEgrav_stars>");
  dump_history_enroll(hst_dEk,"<dEk>");
#ifdef MHD
    dump_history_enroll(hst_dEb,"<dEb>");
    dump_history_enroll(hst_dEb_one,"<dEb_1>");
    dump_history_enroll(hst_B3c,"B3c");
#endif
  dump_history_enroll(hst_Mr_out,"<Mr_out>");
  dump_history_enroll(hst_Mass_out,"<Mass_out>");
  dump_history_enroll(hst_Efree_out,"<Efree_out>");
#ifdef RADIATION
  dump_history_enroll(hst_jsrc,"<jsrc>");
  dump_history_enroll(hst_Er,"<Er>");
  dump_history_enroll(hst_Frho_out,"<Frho_out>");
#endif
#if (NSCALARS > 0)
  dump_history_enroll(hst_MassS0_out,"<MassS0_out>");
  dump_history_enroll(hst_MassS1_out,"<MassS1_out>");
  dump_history_enroll(hst_MassS2_out,"<MassS2_out>");
  dump_history_enroll(hst_MassS3_out,"<MassS3_out>");
  dump_history_enroll(hst_MassS4_out,"<MassS4_out>");
#endif
  dump_history_enroll(hst_mass_gas,"<mass_gas>");
  dump_history_enroll(hst_Reff_out,"\\int R^2 rho dV");
  dump_history_enroll(hst_Meff_out,"\\int rho^2 dV");
  // dump_history_enroll(hst_rhoeffone_out,"\\int rho1 dV");
  // dump_history_enroll(hst_Reffone_out,"\\int R1^2 rho1 dV");
  // dump_history_enroll(hst_Meffone_out,"\\int rho1^2 dV");
  dump_history_enroll(hst_xeff_out,"\\int x rho dV");
  dump_history_enroll(hst_yeff_out,"\\int y rho dV");
  dump_history_enroll(hst_zeff_out,"\\int z rho dV");
  dump_history_enroll(hst_xxeff_out,"\\int xx rho dV");
  dump_history_enroll(hst_yyeff_out,"\\int yy rho dV");
  dump_history_enroll(hst_zzeff_out,"\\int zz rho dV");
  dump_history_enroll(hst_xyeff_out,"\\int xy rho dV");
  dump_history_enroll(hst_xzeff_out,"\\int xz rho dV");
  dump_history_enroll(hst_yzeff_out,"\\int yz rho dV");
  
  /* Enroll BC functions */
  bvals_mhd_fun(pD,left_x1,diode_outflow_ix1);
  bvals_mhd_fun(pD,right_x1,diode_outflow_ox1);
  bvals_mhd_fun(pD,left_x2,diode_outflow_ix2);
  bvals_mhd_fun(pD,right_x2,diode_outflow_ox2);
  bvals_mhd_fun(pD,left_x3,diode_outflow_ix3);
  bvals_mhd_fun(pD,right_x3,diode_outflow_ox3);

#ifdef RADIATION
  /* Set function pointers */
  LuminosityFun = mass_dependent_luminosity;
  kappaP = NULL;
  kappaE = const_absorption;
  kappaF = const_absorption;
#endif
  
  return;
}


/*------------------------------------------------------------------------------
 *  Function ran2
 *
 *  The routine ran2() is extracted from the Numerical Recipes in C
 *  (version 2) code.  I've modified it to use doubles instead of
 *  floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 *  Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 *  with Bays-Durham shuffle and added safeguards.  Returns a uniform
 *  random deviate between 0.0 and 1.0 (exclusive of the endpoint
 *  values).  Call with idum = a negative integer to initialize;
 *  thereafter, do not alter idum between successive deviates in a
 *  sequence.  RNMX should appriximate the largest floating point 
 *  value that is less than 1.
 */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define NTAB 32

double ran2(long int *idum){
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;
  
  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

#undef OFST
#undef KCOMP
#undef KWVM
