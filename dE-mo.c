#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<omp.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

/*********************
*  Global variables  *
*********************/
int fac2(int);
double NormG(double,int,int,int);
double phi(int, int, int, int);
double square(double);

static int Nbasis, Ngaussian, Nshell, Natom, TNatom, HOMOa, HOMOb;
static int *rCrg, *tCrg;
static double *EVa, *EWa, *Gmo, *rCrd, *tCrd;
static double Vr = 0, Vt = 0, dVN = 0;

int main (int argc, char *argv[]){

  int i, j, k;

  if(argc != 3){
    printf("failed, need input arguments\n");
    return;
  }

/*******************************************************
*  read Gaussian data generated from readAO.pl script  *
*******************************************************/
  FILE *fCMNBasis;
  FILE *frCrd, *fCTRC;
  FILE *fEVa, *fEWa, *fACHARGE;
  FILE *fGmo, *fTCrd;

  int tmp;

  /* manipulating input strings */
  char *TCrd;
  char *gmoIn=strdup(argv[1]);

  TCrd=strtok(gmoIn, ".");
  strcat(TCrd,".crd");

  /* read files from perl and python script */
  fCMNBasis = fopen("CMNBasis.pltmp","r");
  frCrd = fopen("COORD.pltmp","r");
  fEVa = fopen("EV_a.pltmp","r");
  fEWa = fopen("EW_a.pltmp","r");
  fGmo = fopen(argv[1],"r");
  fTCrd = fopen(TCrd,"r");
  fACHARGE = fopen("ACHARGE.pltmp","r");


  for(i=0;i<3;i++){
    fscanf(fCMNBasis,"%d",&tmp);
  }
  fscanf(fCMNBasis,"%d",&HOMOa);
  fscanf(fCMNBasis,"%d",&HOMOb);
  fscanf(fCMNBasis,"%d",&Nbasis);
  fscanf(fCMNBasis,"%d",&Ngaussian);
  fscanf(fCMNBasis,"%d",&Nshell);
  fscanf(fCMNBasis,"%d",&Natom);

  /* number of atoms in target system   */
  /* need to be rewrite for alchemy     */
  /* may be done by read until EOF      */
  /* specific format of crd file needed */
  {
    FILE *crdTmp = fopen(TCrd,"r");
    char *string;
    size_t len=0;
    string = (char*) malloc(80*sizeof(char));
    while(!feof(crdTmp)){
      getline(&string,&len,crdTmp);
      TNatom++;
    }
    TNatom--;
  }
  //printf("%d\n",TNatom);

  //printf("%d\t%d\t%d\n",Nbasis,Ngaussian,Nshell);

  EWa = (double *) malloc(Nbasis * sizeof(double));
  EVa = (double *) malloc(Nbasis*Nbasis * sizeof(double));
  Gmo = (double *) malloc(Nbasis*Nbasis * sizeof(double));

  rCrg = (int *) malloc(Natom * sizeof(int));
  rCrd = (double *) malloc(3 * Natom * sizeof(double));
  tCrg = (int *) malloc(TNatom * sizeof(int));
  tCrd = (double *) malloc(3 * TNatom * sizeof(double));

  for(i=0;i<Nbasis;i++){
    fscanf(fEWa,"%lf",&EWa[i]);
    for(j=0;j<Nbasis;j++){
      int s = i*Nbasis + j;
      fscanf(fEVa,"%lf",&EVa[s]);
      fscanf(fGmo,"%lf",&Gmo[s]);
    }
  }

  for(i=0;i<Natom;i++){
    double tmp;
    /* ACHARGE is written in scientific notation */
    fscanf(fACHARGE,"%lf",&tmp);
    rCrg[i] = tmp;
    fscanf(frCrd,"%lf",&rCrd[i*3]);
    fscanf(frCrd,"%lf",&rCrd[i*3+1]);
    fscanf(frCrd,"%lf",&rCrd[i*3+2]);
  }

  for(i=0;i<TNatom;i++){
    fscanf(fTCrd,"%d",&tCrg[i]);
    fscanf(fTCrd,"%lf",&tCrd[i*3]);
    fscanf(fTCrd,"%lf",&tCrd[i*3+1]);
    fscanf(fTCrd,"%lf",&tCrd[i*3+2]);
  }
  //printf("%d %d\n",tCrg[0],rCrg[0]);

  for(i=0;i<Natom;i++){
    for(j=i+1;j<Natom;j++){
      double Rij = 0;
      Rij += square(rCrd[i*3]-rCrd[j*3]);
      Rij += square(rCrd[i*3+1]-rCrd[j*3+1]);
      Rij += square(rCrd[i*3+2]-rCrd[j*3+2]);
      Vr += rCrg[i]*rCrg[j]/sqrt(Rij);
    }
  }
  for(i=0;i<TNatom;i++){
    for(j=i+1;j<TNatom;j++){
      double Rij = 0;
      Rij += square(tCrd[i*3]-tCrd[j*3]);
      Rij += square(tCrd[i*3+1]-tCrd[j*3+1]);
      Rij += square(tCrd[i*3+2]-tCrd[j*3+2]);
      Vt += tCrg[i]*tCrg[j]/sqrt(Rij);
    }
  }
  dVN = Vt - Vr;
  //printf("% 14.7E % 14.7E % 14.7E\n",Vr, Vt, dVN);
/****************END of reading*****************/

  int evi, evj, evk, eva, evb, evc;

  double d1E = 0;
  for(evi=0;evi<HOMOa;evi++){
    int ii = evi*Nbasis + evi;
    d1E += Gmo[ii];
  }
  for(evi=0;evi<HOMOb;evi++){
    int ii = evi*Nbasis + evi;
    d1E += Gmo[ii];
  }
  d1E += dVN;

  double d2E = 0;
  for(evi=0;evi<HOMOa;evi++){
    for(eva=HOMOa;eva<Nbasis;eva++){
      int ia = evi*Nbasis + eva;
      double Cia = EWa[evi] - EWa[eva];
      d2E += square(Gmo[ia]) / Cia;
    }
  }
  for(evi=0;evi<HOMOb;evi++){
    for(eva=HOMOb;eva<Nbasis;eva++){
      int ia = evi*Nbasis + eva;
      double Cia = EWa[evi] - EWa[eva];
      d2E += square(Gmo[ia]) / Cia;
    }
  }

/**************
*  3rd order  *
**************/
  double d3E = 0;
  for(evi=0;evi<HOMOa;evi++){
    for(eva=HOMOa;eva<Nbasis;eva++){
      int ia = evi*Nbasis + eva;
      double Cia = EWa[evi] - EWa[eva];
      for(evb=HOMOa;evb<Nbasis;evb++){
        int ib = evi*Nbasis + evb;
        int ab = eva*Nbasis + evb;
        double Cib = EWa[evi] - EWa[evb];
        double Ciab = Cia*Cib;
        d3E += Gmo[ia]*Gmo[ib]*Gmo[ab] / Ciab;
      }
      for(evj=0;evj<HOMOa;evj++){
        int ja = evj*Nbasis + eva;
        int ij = evj*Nbasis + evj;
        double Cja = EWa[evj] - EWa[eva];
        double Cija = Cia*Cja;
        d3E -= Gmo[ia]*Gmo[ja]*Gmo[ij] / Cija;
      }
    }
  }
  for(evi=0;evi<HOMOb;evi++){
    for(eva=HOMOb;eva<Nbasis;eva++){
      int ia = evi*Nbasis + eva;
      double Cia = EWa[evi] - EWa[eva];
      for(evb=HOMOb;evb<Nbasis;evb++){
        int ib = evi*Nbasis + evb;
        int ab = eva*Nbasis + evb;
        double Cib = EWa[evi] - EWa[evb];
        double Ciab = Cia*Cib;
        d3E += Gmo[ia]*Gmo[ib]*Gmo[ab] / Ciab;
      }
      for(evj=0;evj<HOMOb;evj++){
        int ja = evj*Nbasis + eva;
        int ij = evj*Nbasis + evj;
        double Cja = EWa[evj] - EWa[eva];
        double Cija = Cia*Cja;
        d3E -= Gmo[ia]*Gmo[ja]*Gmo[ij] / Cija;
      }
    }
  }

/**************
*  4th order  *
**************/
  double d4E = 0;

  for(evi=0;evi<HOMOa;evi++){
    for(eva=HOMOa;eva<Nbasis;eva++){
      for(evb=HOMOa;evb<Nbasis;evb++){
        for(evc=HOMOa;evc<Nbasis;evc++){
          int ia = evi*Nbasis + eva;
          int ib = evi*Nbasis + evb;
          int ab = eva*Nbasis + evb;
          int bc = evb*Nbasis + evc;
          double Cia = EWa[evi] - EWa[eva];
          double Cib = EWa[evi] - EWa[evb];
          double Cic = EWa[evi] - EWa[evc];
          double V = Gmo[ia]*Gmo[ib]*Gmo[ab]*Gmo[bc];
          double E = Cia*Cib*Cic;
          d4E += V/E;
        }
      }
    }
  }
 
  for(evi=0;evi<HOMOa;evi++){
    for(evj=0;evj<HOMOa;evj++){
      for(evk=0;evk<HOMOa;evk++){
        for(eva=HOMOa;eva<Nbasis;eva++){
          int ia = evi*Nbasis + eva;
          int ja = evj*Nbasis + eva;
          int ik = evi*Nbasis + evk;
          int jk = evj*Nbasis + evk;
          double Cia = EWa[evi] - EWa[eva];
          double Cja = EWa[evj] - EWa[eva];
          double Cka = EWa[evk] - EWa[eva];
          double V = Gmo[ia]*Gmo[ja]*Gmo[ik]*Gmo[jk];
          double E = Cia*Cja*Cka;
          d4E += V/E;
        }
      }
    }
  }

  for(evi=0;evi<HOMOa;evi++){
    for(evj=0;evj<HOMOa;evj++){
      for(eva=HOMOa;eva<Nbasis;eva++){
        for(evb=HOMOa;evb<Nbasis;evb++){
          int ia = evi*Nbasis + eva;
          int ib = evi*Nbasis + evb;
          int jb = evj*Nbasis + evb;
          int ij = evi*Nbasis + evk;
          int ab = evj*Nbasis + evk;
          double Cia = EWa[evi] - EWa[eva];
          double Cjb = EWa[evj] - EWa[evb];
          double Cib = EWa[evi] - EWa[evb];
          double Cij = EWa[evi] + EWa[evj];
          double Cab = EWa[eva] + EWa[evb];
          double V = Gmo[ia]*Gmo[jb]*Gmo[ij]*Gmo[ab];
          double E = Cia*Cjb*Cib;
          d4E -= 2*V/E;

          double V1 = Gmo[ia]*Gmo[jb]*Gmo[ib]*Gmo[ia];
          double E1 = Cia*Cjb*(Cij - Cab);
          double E2 = Cia*Cjb*Cjb;
          d4E += (V1/E1 -V1/E2);
        }
      }
    }
  }

  d4E *= 2;

//  printf("% 10.7lf\n",d1E);
//  printf("% 10.7lf % 10.7lf\n",d1E, d2E);
//  printf("% 10.7lf % 10.7lf % 10.7f\n",d1E, d2E, d3E);
  printf("% 10.7lf % 10.7lf % 10.7lf % 10.7f\n",d1E, d2E, d3E, d4E);

  free(EVa);
  free(EWa);
  free(Gmo);
  free(tCrg);
  free(rCrg);
  free(tCrd);
  free(rCrd);

  fcloseall();

  return 0;
}

double square(double x){
  return x*x;
}

//int fac2(int n){
//  int val = 1;
//  int k;
//  for(k=n;k>0;k-=2) val *= k;
//  return val;
//}
//
//double NormG(double a, int li, int lj, int lk){
//  int Ni = fac2(2*li-1);
//  int Nj = fac2(2*lj-1);
//  int Nk = fac2(2*lk-1);
//  int Nijk = Ni*Nj*Nk;
//  int nijk = li+lj+lk;
//  double Num1 = pow((2*a/PI),3.0/2.0);
//  double Num2 = pow((4*a),nijk);
//  return pow((Num1*Num2/Nijk),0.5);
//}
//
//double phi(int MOi, int ix, int iy, int iz){
//  int g, s;
//  int evi;
//  int aoi, goi, goI, ai;
//  int px,py,pz;
//  int l;
//  double Rx, Ry, Rz, Ax, Ay, Az, L, L2;
//  double P, a, gR, phi_R = 0;
//  for(evi=0;evi<Nbasis;evi++){
//    /* read indices maps */
//    aoi = AOMAP[evi] - 1;
//    goi = AGMAP[aoi];
//    goI = AGMAP[aoi+1];
//    ai = MAP[aoi] - 1;
//    /* read coordinates and anglular momentum */
//    px = lmMAP[evi*3];
//    py = lmMAP[evi*3+1];
//    pz = lmMAP[evi*3+2];
//    Ax = rCrd[ai*3];
//    Ay = rCrd[ai*3+1];
//    Az = rCrd[ai*3+2];
//    Rx = Sx + ix*dx - Ax;
//    Ry = Sy + iy*dy - Ay;
//    Rz = Sz + iz*dz - Az;
//    L2 = (Rx*Rx + Ry*Ry + Rz*Rz);
//    P = pow(Rx,px) * pow(Ry,py) * pow(Rz,pz);
//    s = MOi*Nbasis + evi;
//    for(g=goi;g<goI;g++){
//      /* evaluate gaussian functions */
//      a = PEXP[g];
//      gR = CTRC[g]*NormG(a,px,py,pz)*P*exp(-a*L2)*EVa[s];
//      if(fabs(gR) > pow(10,-15)) phi_R += gR;
//    }
//  }
//  return phi_R;
//}
