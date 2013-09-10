#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <iostream>
//#include <windows.h>
#include<fstream>
#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include <time.h> 
#include <cstdlib>   // for rand
#include <cmath>     // for atan, sqrt, log, cos
#include <algorithm> // for generate_n
using namespace std;


#define NDIM 3
//#define MaxInt 21474836

typedef struct {double x, y, z;} VecR;
typedef struct {int x, y, z;} VecI;



#define Sqr(x)     ((x) * (x))
#define Cube(x)    ((x) * (x) * (x))
#define Sgn(x, y)  (((y) >= 0) ? (x) : (- (x)))
#define IsEven(x)  ((x) & ~1)
#define IsOdd(x)   ((x) & 1)
#define Nint(x)    (((x)<0.) ? (-(int)(0.5- (x))): ((int)(0.5 + (x))))
#define Min(x1, x2)  (((x1) < (x2)) ? (x1) : (x2))
#define Max(x1, x2)  (((x1) > (x2)) ? (x1) : (x2))
#define Min3(x1, x2, x3) \
   (((x1) < (x2)) ? (((x1) < (x3)) ? (x1) : (x3)) :         \
                    (((x2) < (x3)) ? (x2) : (x3)))
#define Max3(x1, x2, x3) \
   (((x1) > (x2)) ? (((x1) > (x3)) ? (x1) : (x3)) :         \
                    (((x2) > (x3)) ? (x2) : (x3)))
#define Clamp(x, lo, hi)                                    \
   (((x) >= (lo) && (x) <= (hi)) ? (x) :                    \
   (((x) < (lo)) ? (lo) : (hi)))




#define VSet(v, sx, sy, sz)                                 \
   (v).x = sx,                                              \
   (v).y = sy,                                              \
   (v).z = sz
#define VCopy(v1, v2)                                       \
   (v1).x = (v2).x,                                         \
   (v1).y = (v2).y,                                         \
   (v1).z = (v2).z
#define VScale(v, s)                                        \
   (v).x *= s,                                              \
   (v).y *= s,                                              \
   (v).z *= s
#define VSCopy(v2, s1, v1)                                  \
   (v2).x = (s1) * (v1).x,                                  \
   (v2).y = (s1) * (v1).y,                                  \
   (v2).z = (s1) * (v1).z
#define VAdd(v1, v2, v3)                                    \
   (v1).x = (v2).x + (v3).x,                                \
   (v1).y = (v2).y + (v3).y,                                \
   (v1).z = (v2).z + (v3).z
#define VSub(v1, v2, v3)                                    \
   (v1).x = (v2).x - (v3).x,                                \
   (v1).y = (v2).y - (v3).y,                                \
   (v1).z = (v2).z - (v3).z
#define VMul(v1, v2, v3)                                    \
   (v1).x = (v2).x * (v3).x,                                \
   (v1).y = (v2).y * (v3).y,                                \
   (v1).z = (v2).z * (v3).z
#define VDiv(v1, v2, v3)                                    \
   (v1).x = (v2).x / (v3).x,                                \
   (v1).y = (v2).y / (v3).y,                                \
   (v1).z = (v2).z / (v3).z
#define VSAdd(v1, v2, s3, v3)                               \
   (v1).x = (v2).x + (s3) * (v3).x,                         \
   (v1).y = (v2).y + (s3) * (v3).y,                         \
   (v1).z = (v2).z + (s3) * (v3).z
#define VSSAdd(v1, s2, v2, s3, v3)                          \
   (v1).x = (s2) * (v2).x + (s3) * (v3).x,                  \
   (v1).y = (s2) * (v2).y + (s3) * (v3).y,                  \
   (v1).z = (s2) * (v2).z + (s3) * (v3).z
#define VDot(v1, v2)                                        \
   ((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)
#define VWDot(v1, v2, v3)                                   \
   ((v1).x * (v2).x * (v3).x + (v1).y * (v2).y * (v3).y +   \
   (v1).z * (v2).z * (v3).z)
#define VCross(v1, v2, v3)                                  \
   (v1).x = (v2).y * (v3).z - (v2).z * (v3).y,              \
   (v1).y = (v2).z * (v3).x - (v2).x * (v3).z,              \
   (v1).z = (v2).x * (v3).y - (v2).y * (v3).x
#define MVMul(v1, m, v2)                                         \
   (v1).x = (m)[0] * (v2).x + (m)[3] * (v2).y + (m)[6] * (v2).z, \
   (v1).y = (m)[1] * (v2).x + (m)[4] * (v2).y + (m)[7] * (v2).z, \
   (v1).z = (m)[2] * (v2).x + (m)[5] * (v2).y + (m)[8] * (v2).z
#define MVMulT(v1, m, v2)                                        \
   (v1).x = (m)[0] * (v2).x + (m)[1] * (v2).y + (m)[2] * (v2).z, \
   (v1).y = (m)[3] * (v2).x + (m)[4] * (v2).y + (m)[5] * (v2).z, \
   (v1).z = (m)[6] * (v2).x + (m)[7] * (v2).y + (m)[8] * (v2).z
#define VProd(v)                                            \
   ((v).x * (v).y * (v).z)
#define VGe(v1, v2)                                         \
   ((v1).x >= (v2).x && (v1).y >= (v2).y && (v1).z >= (v2).z)
#define VLt(v1, v2)                                         \
   ((v1).x < (v2).x && (v1).y < (v2).y && (v1).z < (v2).z)
#define VLinear(p, s)                                       \
   (((p).z * (s).y + (p).y) * (s).x + (p).x)
#define VSetAll(v, s)                                       \
   VSet (v, s, s, s)
#define VAddCon(v1, v2, s)                                  \
   (v1).x = (v2).x + (s),                                   \
   (v1).y = (v2).y + (s),                                   \
   (v1).z = (v2).z + (s)
#define VComp(v, k)                                         \
   *((k == 0) ? &(v).x : ((k == 1) ? &(v).y : &(v).z))
#define VToLin(a, n, v)                                     \
   a[(n) + 0] = (v).x,                                      \
   a[(n) + 1] = (v).y,                                      \
   a[(n) + 2] = (v).z
#define VFromLin(v, a, n)                                   \
   VSet (v, a[(n) + 0], a[(n) + 1], a[(n) + 2])
#define VCSum(v)                                            \
   ((v).x + (v).y + (v).z)



#define VZero(v)  VSetAll (v, 0)
#define VLenSq(v)  VDot (v, v)
#define VWLenSq(v1, v2)  VWDot(v1, v2, v2)
#define VLen(v)  sqrt (VDot (v, v))
#define VVAdd(v1, v2)  VAdd (v1, v1, v2)
#define VVSub(v1, v2)  VSub (v1, v1, v2)
#define VVSAdd(v1, s2, v2) VSAdd (v1, v1, s2, v2)
#define VInterp(v1, s2, v2, v3)                             \
   VSSAdd (v1, s2, v2, 1. - (s2), v3)

typedef struct {
   double u1, u2, u3, u4;
} Quat;

#define QSet(q, s1, s2, s3, s4)                             \
   (q).u1 = s1,                                             \
   (q).u2 = s2,                                             \
   (q).u3 = s3,                                             \
   (q).u4 = s4
#define QZero(q)  QSet (q, 0, 0, 0, 0)
#define QScale(q, s)                                        \
   (q).u1 *= s,                                             \
   (q).u2 *= s,                                             \
   (q).u3 *= s,                                             \
   (q).u4 *= s
#define QSAdd(q1, q2, s3, q3)                               \
   (q1).u1 = (q2).u1 + (s3) * (q3).u1,                      \
   (q1).u2 = (q2).u2 + (s3) * (q3).u2,                      \
   (q1).u3 = (q2).u3 + (s3) * (q3).u3,                      \
   (q1).u4 = (q2).u4 + (s3) * (q3).u4
#define QLenSq(q)                                           \
   (Sqr ((q).u1) + Sqr ((q).u2) + Sqr ((q).u3) +            \
   Sqr ((q).u4))
#define QMul(q1, q2, q3)                                    \
   (q1).u1 =   (q2).u4 * (q3).u1 - (q2).u3 * (q3).u2 +      \
               (q2).u2 * (q3).u3 + (q2).u1 * (q3).u4,       \
   (q1).u2 =   (q2).u3 * (q3).u1 + (q2).u4 * (q3).u2 -      \
               (q2).u1 * (q3).u3 + (q2).u2 * (q3).u4,       \
   (q1).u3 = - (q2).u2 * (q3).u1 + (q2).u1 * (q3).u2 +      \
               (q2).u4 * (q3).u3 + (q2).u3 * (q3).u4,       \
   (q1).u4 = - (q2).u1 * (q3).u1 - (q2).u2 * (q3).u2 -      \
               (q2).u3 * (q3).u3 + (q2).u4 * (q3).u4

typedef struct {
  double R, I;
} Cmplx;

#define CSet(a, x, y)                                       \
   a.R = x,                                                 \
   a.I = y
#define CAdd(a, b, c)                                       \
   a.R = b.R + c.R,                                         \
   a.I = b.I + c.I
#define CSub(a, b, c)                                       \
   a.R = b.R - c.R,                                         \
   a.I = b.I - c.I
#define CMul(a, b, c)                                       \
  a.R = b.R * c.R - b.I * c.I,                              \
  a.I = b.R * c.I + b.I * c.R




/********** End of All definitive defines**********#include "vdefs.h"*/






//#include "namelist.h"

typedef enum {N_I, N_R} VType;

#define NameI(x)  {#x, &x, N_I, sizeof (x) / sizeof (int)}
#define NameR(x)  {#x, &x, N_R, sizeof (x) / sizeof (double)}

typedef struct {
  char *vName;
  void *vPtr;
  VType vType;
  int vLen, vStatus;
} NameList;

#define ValI(x)  {&x, N_I, sizeof (x) / sizeof (int)}
#define ValR(x)  {&x, N_R, sizeof (x) / sizeof (double)}

typedef struct {
  void *vPtr;
  VType vType;
  int vLen;
} ValList;




//typedef double double;

#define CHAR_MINUS  '-'
#define CHAR_ZERO   '0'
#define M_PI 3.14159265358979323846
//#define NDIM  3


typedef struct {
  double u[9];
} RMat;

#define MAT(a, n, i, j)  (a)[(i) + n * (j)]

#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))

#define AllocMem2(a, n1, n2, t)                             \
   AllocMem (a, n1, t *);                                   \
   AllocMem (a[0], (n1) * (n2), t);                         \
   for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;

#define MAX_MPEX_ORD  2
#define I(i, j)  ((i) * ((i) + 1) / 2 + (j))
#define c(i, j)  c[I(i, j)]
#define s(i, j)  s[I(i, j)]



#define DO_MOL  for (n = 0; n < nMol; n ++)
#define DO_CELL(j, m)  for (j = cellList[m]; j >= 0; j = cellList[j])


/*
#define VWrap(v, t)                                         \
   if (v.t >= 0.5 * region.t)      v.t -= region.t;         \
   else if (v.t < -0.5 * region.t) v.t += region.t

#define VShift(v, t)                                        \
   if (v.t >= 0.5 * region.t)      shift.t -= region.t;     \
   else if (v.t < -0.5 * region.t) shift.t += region.t



#define VShiftWrap(v, t)                                    \
   if (v.t >= 0.5 * region.t) {                             \
     shift.t -= region.t;                                   \
     v.t -= region.t;                                       \
   } else if (v.t < -0.5 * region.t) {                      \
     shift.t += region.t;                                   \
     v.t += region.t;                                       \
   }


#define VCellWrap(t)                                        \
   if (m2v.t >= cells.t) {                                  \
     m2v.t = 0;                                             \
     shift.t = region.t;                                    \
   } else if (m2v.t < 0) {                                  \
     m2v.t = cells.t - 1;                                   \
     shift.t = - region.t;                                  \
   }




*/



#define VWrap(v, t)                                      \
   if (v.t >= 0.5 * simulationregion.t)      v.t -= simulationregion.t;         \
   else if (v.t < -0.5 * simulationregion.t) v.t += simulationregion.t



#define VShift(v, t)                                        \
   if (v.t >= 0.5 * simulationregion.t)      shift.t -= simulationregion.t;     \
   else if (v.t < -0.5 * simulationregion.t) shift.t += simulationregion.t


#define VShiftWrap(v, t)                                    \
   if (v.t >= 0.5 * simulationregion.t) {                             \
     shift.t -= simulationregion.t;                                   \
     v.t -= simulationregion.t;                                       \
   } else if (v.t < -0.5 * simulationregion.t) {                      \
     shift.t += simulationregion.t;                                   \
     v.t += simulationregion.t;                                       \
   }


#define VCellWrap(t)                                        \
   if (m2v.t >= cells.t) {                                  \
     m2v.t = 0.;                                             \
     shift.t =simulationregion.t;                                    \
   } else if (m2v.t < 0) {                                  \
     m2v.t = cells.t - 1.;                                   \
     shift.t = - simulationregion.t;                                  \
   }

/*adding by wu*/

#define VWrapR (v) v-2*v/VLenSq(v)




#define VWrapAll(v)                                         \
   {VWrap (v, x);                                           \
   VWrap (v, y);                                            \
   VWrap (v, z);}
#define VShiftAll(v)                                        \
   {VShift (v, x);                                          \
   VShift (v, y);                                           \
   VShift (v, z);}
#define VCellWrapAll()                                      \
   {VCellWrap (x);                                          \
   VCellWrap (y);                                           \
   VCellWrap (z);}


#define OFFSET_VALS                                           \
   { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},            \
     {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},  \
     {-1,-1,1}, {0,-1,1}, {1,-1,1}                            \
   }

#define N_OFFSET  14





#define HalfVCellWrap(t)                                        \
   if (m2v.t >= halfdrop_cells.t) {                                  \
     m2v.t = 0;                                             \
     shift.t =simulationregion.t;                                    \
   } else if (m2v.t < 0) {                                  \
     m2v.t = halfdrop_cells.t - 1.;                                   \
     shift.t = - simulationregion.t;                                  \
   }


typedef struct {
  double val, sum, sum2;
} Prop;

#define PropZero(v)  v.sum = v.sum2 = 0.
#define PropAccum(v)  v.sum += v.val, v.sum2 += Sqr (v.val)
#define PropAvg(v, n) \
   v.sum /= n, v.sum2 = sqrt (Max (v.sum2 / n - Sqr (v.sum), 0.))
#define PropEst(v)  v.sum, v.sum2

typedef struct {
  double time;
  int left, right, up, circAL, circAR, circBL, circBR, idA, idB;
} EvTree;

#define MOL_LIMIT  10000000

#define NameVal(x)                                          \
   if (! strncmp (bp, #x, strlen (#x))) {                   \
     bp += strlen (#x);                                     \
     x = strtod (bp, &bp);                                  \
   }





#define NHIST  (NDIM + 2)

#define ReadF(x)       fread  (&x, sizeof (x), 1, fp)
#define WriteF(x)      fwrite (&x, sizeof (x), 1, fp)
#define ReadFN(x, n)   fread  (x, sizeof (x[0]), n, fp)
#define WriteFN(x, n)  fwrite (x, sizeof (x[0]), n, fp)

enum {FL_CHECKA, FL_CHECKB, FL_CKLAST, FL_SNAP};
char *fileNameR[] = {"xxnnchecka.data", "xxnncheckb.data",
   "xxnncklast.data", "xxnnsnap.data"}, fileName[5][20];

char *progId = "md";

enum {ERR_NONE, ERR_BOND_SNAPPED, ERR_CHECKPT_READ, ERR_CHECKPT_WRITE,
   ERR_COPY_BUFF_FULL, ERR_EMPTY_EVPOOL, ERR_MSG_BUFF_FULL,
   ERR_OUTSIDE_REGION, ERR_SNAP_READ, ERR_SNAP_WRITE,
   ERR_SUBDIV_UNFIN, ERR_TOO_MANY_CELLS, ERR_TOO_MANY_COPIES,
   ERR_TOO_MANY_LAYERS, ERR_TOO_MANY_LEVELS, ERR_TOO_MANY_MOLS,
   ERR_TOO_MANY_MOVES, ERR_TOO_MANY_NEBRS, ERR_TOO_MANY_REPLICAS};

char *errorMsg[] = {"", "bond snapped", "read checkpoint data",
   "write checkpoint data", "copy buffer full", "empty event pool",
   "message buffer full", "outside region", "read snap data",
   "write snap data", "subdivision unfinished", "too many cells",
   "too many copied mols", "too many layers", "too many levels",
   "too many mols", "too many moved mols", "too many neighbors",
   "too many replicas"};


typedef struct {
  VecR r, rv, ra, ra1, ra2, ro, rvo;
} Mol;

Mol *mol;
Mol *halfmol;//adding by wu



VecR region, vSum,fSum;
VecR simulationregion;
VecI initUcell={4,40,40};
double deltaT=0.002, density=0.4, rCut=3.8, temperature=0.70758, timeNow, uSum, velMag, vvSum,virSum;
double regionR;//adding by wu
Prop kinEnergy, totEnergy,pressure;
int moreCycles, nMol, stepAvg=2000, stepCount, stepEquil=1000, stepLimit=2500000, runId=1,stepOutPut=100;

VecI cells;//cells (x,y,z)::the number of cell division of physical region=cells.x *cells.y*cells.z;
           //its value is assigned in SetParams()
int *cellList;//linked list 

double dispHi, rNebrShell=0.4;
int *nebrTab, nebrNow, nebrTabFac=8, nebrTabLen, nebrTabMax;
VecR obsPos={-0.25,0.0};

VecI sizeHistGrid={1,1,50};
double **histGrid, bdyStripWidth=3.0, flowSpeed=2.0, obsSize=0.2;
double gravField=0.1;
int countGrid, limitGrid=100, nFixedMol, nFreeMol, snapNumber, stepDrive=40,
   stepGrid=500;

int HalfNmol;//adding by Wu
double *profileT,*profileV;


//velocity distribution//

double *histVel, rangeVel=0.0;
int countVel, limitVel=4,sizeHistVel=50, stepVel=5;


//end of velocity distribution




/*New added parameter inputs by user*/
//double  bdyStripWidth=3.0, flowSpeed=2.0, obsSize=0.2;


/*
bdyStripWidth=3.0;
deltaT=0.005;
density=0.8;
flowSpeed=2.0;
iniUcell={500,250};
limitGrid=200;
nebrTabFac=8;
obsPos={-0.25,0.0}
obsSize=0.2;
rNebrShell=0.4;
runId=1;
sizeHistGrid={120,60};
stepAvg=100;
stepDrive=40;
stepEquil=0;
stepGrid=20;
stepLimit=100000;
temperature=1.0;
*/
/********************88888888888*/

NameList nameList[] = {
  NameR (bdyStripWidth),
  NameR (deltaT),
  NameR (density),
  NameR (flowSpeed),
  NameI (initUcell),
  NameI (limitGrid),
  NameI (nebrTabFac),
  NameR (obsPos),
  NameR (obsSize),
  NameR (rNebrShell),
  NameI (runId),
  NameI (sizeHistGrid),
  NameI (stepAvg),
  NameI (stepDrive),
  NameI (stepEquil),
  NameI (stepGrid),
  NameI (stepLimit),
  NameR (temperature),
};



/***************To generate random number and vector**************/





/*
void InitRand (int randSeedI)//based on the system clock to initialize the random number sequence
{
  
   SYSTEMTIME tv;

   if(randSeedI != 0) randSeedP = randSeedI;
   else {
    GetLocalTime(&tv);
    randSeedP = tv.wSecond;
  }
}*/

 int randSeedP=25;


double RandR ()
{
	
#define SCALE_FAC  32767.

#define IADD   283806245
#define IMUL   314159269
#define MASK   2147483647
#define SCALE  0.4656612873e-9

//int randSeedP= rand() % 10 + stepCount;

//int randSeedP = 100;//for reproducible run 
  randSeedP = (randSeedP * IMUL + IADD) & MASK;
  
  if ((randSeedP * SCALE)>1)
  {
	  exit(0);
  }
  return (randSeedP * SCALE);
   
}

#if NDIM == 2

void VRand (VecR *p)
{
  double s;

  s = 2. * M_PI * RandR ();
  p->x = cos (s);
  p->y = sin (s);
}

#elif NDIM == 3
/*
void VRand (VecR *p)
{
  double s, x, y;

  s = 2.;
  while (s > 1.) {
    x = 2. * RandR () - 1.;
    y = 2. * RandR () - 1.;
    s = Sqr (x) + Sqr (y);
  }
  p->z = 1. - 2. * s;
  s = 2. * sqrt (1. - s);
  p->x = s * x;
  p->y = s * y;
}
*/

double gasdev(double sigma, double mue) 
{
static bool available = false;
static double gset;
double fac, rsq, v1, v2;
if (!available) {
do {
//srand( time( NULL ) );
v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
//srand ( time(NULL) );
v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
rsq = v1 * v1 + v2 * v2;

} while (rsq >= 1.0 || rsq == 0.0);
fac = sqrt(-2.0 * log(rsq) / rsq);
gset = v1 * fac;
 
available = true;
return v2*fac;
} else {
available = false;
gset=gset*sigma+mue;
return gset;
}
}


void VRand (VecR *p,double sigma)
{
  double s, x, y;

  s = 2.;
  while (s > 1.) {
    x = 2. * gasdev (sigma,0) - 1.;
    y = 2. * gasdev (sigma,0) - 1.;
	 //   x = 2. * RandR () - 1.;
  //  y = 2. * RandR () - 1.;

    s = Sqr (x) + Sqr (y);
  }
  p->z = 1. - 2. * s;
  s = 2. * sqrt (1. - s);
  p->x = s * x;
  p->y = s * y;
}


#endif









void ErrExit (int code)
{
  printf ("Error: %s\n", errorMsg[code]);
  exit (0);
}

/***************To generate random number and vector**************/

//to include all declarations for functions#include "prototype.h"

void AccumBondAngDistn (int);
void AccumDiffusion (void);
void AccumDihedAngDistn (int);
void AdjustInitTemp (void);
void AccumProps (int);
void AccumSpacetimeCorr (void);
void AccumVacf (void);
void AddBondedPair (int, int);
void AdjustDipole (void);
void AdjustLinkAngles (void);
void AdjustPressure (void);
void AdjustQuat (void);
void AdjustTemp(void);

int stepAdjustTemp=1000;
void AdjustVelocity (void);
void AllocArrays (void);
void AnalClusterSize (void);
void AnalVorPoly (void);
void AnlzConstraintDevs (void);
void ApplyBarostat (void);
void ApplyBoundaryCond (void);
void *ApplyBoundaryCondT (void *);
void ApplyThermostat (void);
void ApplyWallBoundaryCond (void);
void AssignMpCells (void);
void AssignToChain (void);
void BisectPlane (void);
void BuildClusters (void);
void BuildConstraintMatrix (void);
void BuildIntTree (void);
void BuildLinkInertiaMats (void);
void BuildLinkMmat (double *, int);
void BuildLinkPhimatT (double *, int);
void BuildLinkRotmatT (RMat *, double, double);
void BuildLinkXYvecs (int);
void BuildNebrList (void);
void *BuildNebrListT (void *);
void BuildRotMatrix (RMat *, Quat *, int);
void BuildStepRmatT (RMat *, VecR *);
void CombineMpCell (void);
void CompressClusters (void);
void ComputeAccelsQ (void);
void ComputeAngVel (int, VecR *);
void ComputeBdyForces (void);
void ComputeChainAngleForces (void);
void ComputeChainBondForces (void);
void ComputeChainTorsionForces (void);
void ComputeConstraints (void);
void ComputeDerivsPT (void);
void ComputeDipoleAccel (void);
void ComputeExternalForce (void);
void ComputeFarCellInt (void);


void ComputeForces (void);
void ComputeForcesNebr (void);
void ComputeForcesCellDiv(void);

void CoolingComputeForcesCellDiv(void);
void relaxApplyThermostat (void);
void ComputeForcesPair(void);


void ComputeForcesDipoleF (void);
void ComputeForcesDipoleR (void);
void *ComputeForcesT (void *);
void ComputeLinkCoordsVels (void);
void ComputeLinkAccels (void);
void ComputeLinkForces (void);
void ComputeNearCellInt (void);
void ComputeSiteForces (void);
void ComputeThermalForce (void);
void ComputeTorqs (void);
void ComputeWallForces (void);
void CorrectorStep (void);
void CorrectorStepBox (void);
void CorrectorStepPT (void);
void CorrectorStepQ (void);
void CorrectorStepS (void);
void CrystalInitCoords (void);
void DefineMol (void);
void DeleteAllMolEvents (int id);
void DeleteEvent (int);
void DoPackInt (int *, int);
void DoPackReal (double *, int);
void DoParlCopy (void);
void DoParlMove (void);
void DoUnpackInt (int *, int);
void DoUnpackReal (double *, int);
void DriveFlow (void);
void ErrExit (int);
void EulerToQuat (Quat *, double *);
void EvalMpCell (void);
void EvalMolCount (void);
void EvalChainProps (void);
void EvalDiffusion (void);
void EvalDihedAngCorr (void);
void EvalEamParams (void);
void EvalFreePath (void);
void EvalHelixOrder (void);
void EvalLatticeCorr (void);
//void EvalMpForce (VecR *, double *, MpTerms *, MpTerms *, int);
//void EvalMpL (MpTerms *, VecR *, int);
//void EvalMpM (MpTerms *, VecR *, int);
//void EvalMpProdLL (MpTerms *, MpTerms *, MpTerms *, int);
//void EvalMpProdLM (MpTerms *, MpTerms *, MpTerms *, int);
void EvalProfile (void);
void EvalProps (void);
void EvalRdf (void);
void EvalSinCos (void);
void EvalSpacetimeCorr (void);
void EvalVacf (void);
void EvalVelDist (void);
void FftComplex (Cmplx *, int);
void FindDistVerts (void);
void FindTestSites (int);
void GatherWellSepLo (void);
void GenSiteCoords (void);
void GetCheckpoint (void);
int  GetConfig (void);
int  GetGridAverage (void);

void  GetNameList ();
void GridAverage (int);
void InitAccels (void);
void InitAngAccels (void);
void InitAngCoords (void);
void InitAngVels (void);
void InitBoxVars (void);
void InitCharges (void);
void InitClusters (void);
void CookingCoords (void);
void PrintInitCoord(void);//output of initcoordinate adding by wu
void PrintCookingCoord(void);//output of initcoordinate adding by wu
void PrintCoordInit(void);//output of latercoordinate adding by wu

void coolingprocess(void);
void CoolingStep(void);
void relaxartificialforce(void);
void relaxStep ( int );
void ComputeForcesCellDiv3(int);
void PrintDropArtificialCoord(void);
void PrintDropCoord(void);


void HemisphereInit(void);// adding by wu
void PrinthalfdropletCoordInit(void);//adding by wu

void PrintVeloInit(void);//output of initial velocity by wu



void InitCoordsWalls (double);
void InitDiffusion (void);
void InitEventList (void);
void InitFeedbackVars (void);
void InitFreePath (void);
void InitLinkState (void);
void InitPairEng (void);
void InitRand (int);
void InitSlaves (void);
void InitSpacetimeCorr (void);
void InitState (void);
void InitVacf (void);
void InitVels (void);
void InitVorPoly (void);
double Integrate (double*, int);
void LeapfrogStep (int);
void LeapfrogStepLinks (int);
void *LeapfrogStepT (void *);
void LocateIntTreeCellCm (void);
void MeasureTrajDev (void);
void MulMat (double *, double *, double *, int);
void MulMatVec (double *, double *, double *, int);
void MultipoleCalc (void);
void NebrParlProcs (void);
void NextEvent (void);

int OutsideObs (VecR *);
void PackCopiedData (int, int, int *, int);
void PackMovedData (int, int, int *, int);
void PackValList (ValList *, int);
void PerturbCoords (void);
void PerturbTrajDev (void);
void PolyGeometry (void);
void PolySize (void);
void PredictEvent (int, int);
void PredictorStep (void);
void PredictorStepBox (void);
void PredictorStepPT (void);
void PredictorStepQ (void);
void PredictorStepS (void);
void PrintChainProps (FILE *);
void PrintDiffusion (FILE *);
void PrintDihedAngCorr (FILE *);
void PrintFreePath (FILE *);
void PrintHelp (char *);
void PrintNameList ();
void PrintPairEng (FILE *);
void PrintProfile ();
void PrintRdf (FILE *);
void PrintSpacetimeCorr (FILE *);
void PrintSummary ();
void PrintTrajDev (FILE *);
void PrintVacf (FILE *);
void PrintVelDist (void);
void ProcCutEdges (void);
void ProcCutFaces (void);
void ProcDelVerts (void);
void ProcessCellCrossing (void);
void ProcessCollision (void);
void ProcInterrupt ();
void PropagateCellLo (void);
void ProcNewFace (void);
void ProcNewVerts (void);
void PutCheckpoint (void);
void PutConfig (void);
void PutPlotData (void);
void PutGridAverage (void);
double RandR (void);
void RemoveOld (void);
void RepackMolArray (void);
void ReplicateMols (void);
void RestoreConstraints (void);
void ScaleCoords (void);
void ScaleVels (void);
void ScanIntTree (void);
void ScheduleEvent (int, int, double);
void SetMolType (void);
void SetBase (void);
void SetCellSize (void);
void SetMolSizes (void);
void SetParams (void);
void SetupFiles (void);
void SetupInterrupt (void);
void SetupJob (void);
void SetupLayers (void);
void SingleEvent (void);
void SingleStep (void);
void SolveCubic (double *, double *);
void SolveLineq (double *, double *, int);
void Sort (double *, int *, int);
void StartRun (void);
void SubdivCells (void);
void UnpackCopiedData (int);
void UnpackMovedData (int);
void UnpackValList (ValList *, int);
void UnscaleCoords (void);
void UpdateMol (int);
void UpdateCellSize (void);
void UpdateSystem (void);
void VRand (VecR *, double);
void ZeroDiffusion (void);
void ZeroFixedAccels (void);
void ZeroSpacetimeCorr (void);
void ZeroVacf (void);


/************End of decrealation********************/

/**********decrealation associated with solid substrate****************/


#define DO_Solid_MOL  for (n = 0; n < solid_nMol; n ++)
#define DO_Solid_CELL(j, m)  for (j = solid_cellList[m]; j >= 0; j = solid_cellList[j])


void SetParams_solid(void);

double solid_density;
int solid_nMol;
int *solid_cellList;//linked list
Mol *solid_mol;

VecR solid_region;
VecR solid_simulationregion, solid_initUcell={12,12,2};
VecI solid_cells;
void Solid_SetupJob (void);



void CrystalSolidCoords(void);
void PrintSolidCoord(void);//initCoords output
void SolidInitAccels (void);
void SolidInitVels (void);
void SolidSingleStep(void);

void SolidLeapfrogStep (int);
void SolidApplyBoundaryCond(void);
void SolidComputeForcesCellDiv(void);
void PrintSolidCoordInit(void);
  	

#define SolidVCellWrap(t)                                        \
   if (m2v.t >= solid_cells.t) {                                  \
     m2v.t = 0;                                             \
     shift.t = simulationregion.t;                                    \
   } else if (m2v.t < 0) {                                  \
     m2v.t = solid_cells.t - 1;                                   \
     shift.t = - simulationregion.t;                                  \
   }


#define SolidVCellWrapAll()                                      \
   {SolidVCellWrap (x);                                          \
   SolidVCellWrap (y);                                           \
   SolidVCellWrap (z);}


//assemble system
int system_nMol;
#define DO_System_MOL  for (n = 0; n < system_nMol; n ++)
#define DO_System_CELL(j, m)  for (j = system_cellList[m]; j >= 0; j = system_cellList[j])

Mol *system_mol;
int *system_cellList;//linked list

void solid_drop_systemInit(void);
void PrintSystemCoordInit(void);
void PrintEquilSystemCoord(void);//initial system
void SystemMove(void);

//to study system of drop spreading

void System_SingleStep (void);
void System_LeapfrogStep (int);
void System_ApplyBoundaryCond(void);
void System_ComputeForcesCellDiv(void);
VecI system_cells;
VecR simulationsystemregion;
double topwall_z=0.0;
 

#define VWrapAllXY(v)                                         \
   {VWrap (v, x);                                           \
   VWrap (v, y);                                            \
   }

#define VWrapAllZ(v)                                         \
   {VWrapZ (v, z);                                           \
                                     \
   }


//*************************structureless simulation******************/

VecI halfdrop_cells;
int *halfdrop_cellList;//linked list
void dropspreading(void);
void Structureless_SystemMove(void);
void Structureless_SingleStep (void);
void PrintEquilStructurelessCoord(void);

void Structureless_LeapfrogStep (int);
void Structureless_ApplyBoundaryCond(void);//(the periodic conditions work for x,y directions; but not z)working on here
void Structureless_ComputeForcesCellDiv(void);//working on here
double NHamaker=0.2;
void HalfdropInitVels (void);
void HalfdropInitAccels(void);//ok
void PrintStructurelessCoordEvalution(void);
void PrintStructurelessCoordEvalution2(void);

void CoolingApplyThermostat (void);
void extremalzvalue(void);

//************************End of structureless simulation******************/


void DropInitReading(void);
void PrintEquilStructurelessAcc(void);
void Structureless_ComputeForcesPair(void);
void PrintEquilCoordVelAcc(void);
void Printdropspreadingprocess(void);

#define PCR4(r, ro, v, a, a1, a2, t)                        \
   r.t = ro.t + deltaT * v.t +                              \
   wr * (cr[0] * a.t + cr[1] * a1.t + cr[2] * a2.t)
#define PCV4(r, ro, v, a, a1, a2, t)                        \
   v.t = (r.t - ro.t) / deltaT +                            \
   wv * (cv[0] * a.t + cv[1] * a1.t + cv[2] * a2.t)


#define PR(t)                                               \
   PCR4 (halfmol[n].r, halfmol[n].r, halfmol[n].rv,                     \
   halfmol[n].ra, halfmol[n].ra1, halfmol[n].ra2, t)
#define PRV(t)                                              \
   PCV4 (halfmol[n].r, halfmol[n].ro, halfmol[n].rv,                    \
   halfmol[n].ra, halfmol[n].ra1, halfmol[n].ra2, t)
#define CR(t)                                               \
   PCR4 (halfmol[n].r,halfmol[n].ro, halfmol[n].rvo,                   \
   halfmol[n].ra, halfmol[n].ra1, halfmol[n].ra2, t)
#define CRV(t)                                              \
   PCV4 (halfmol[n].r, halfmol[n].ro, halfmol[n].rv,                    \
   halfmol[n].ra, halfmol[n].ra1, halfmol[n].ra2, t)





// Pattern substrate
double Radius=initUcell.y;
//double HWGR=0.3;//the ratio of H/R, W/R, G/R
//double Multiplier=0.5;


/*int NLAYERS=int(HWGR*Radius);

double StepGap=HWGR*Radius;
double StepWidth=HWGR*Radius;
int NSTEPS=int((6*Radius)/(StepGap+StepWidth));*/


double StepWidth=0.2*Radius;
double StepHeight=0.25*Radius;
double StepGap=0.45*Radius;
int    NSTEPS=int((7*Radius+StepGap)/(2*(StepGap+StepWidth)))*2;
double Lgs=0.921; //(1+0.3)/2 Lgs=0.95 from Bojan //0.3--penetrating into sidewall
double Egs=0.25; //Change Energy to vary wetting behavior
double rous=4.0/powl(1.2,3);
double Las=sqrtl(3.0); //area  of lattice 




double InitialMove=0.8;//initial distance above the top step

double duattractive (double , double );
double durepulsive (double , double );
double duattractiveapp (double , double );
double durepulsiveapp (double , double );
double Sign(double );

void PrinthalfdropletCoordMove(void);
double durepulsiveyLz (double , double );//y>>z
double duattractiveyLz (double , double );
double durepulsiveySz (double , double );//y<<z
double duattractiveySz (double , double );

double durepulsiveyEz (double , double );//y~O(z)
double duattractiveyEz (double , double );



void PrintStructurelessForceEvalution(void);
double sidewallforce(int, double ,double ,double );
void forcenormalization(void);
double fcValmaxw=0.0;
double fcValmaxd=0.0;
double fcValmaxupperwall=0.0;
double fcValmaxstep=0.0;
void ApplyThermostat2(void);
double cdr=0.4;
double cda=1.0;
double CLA=1.0;
double CLL=1.0;
double epsilon=1.0e-9;
void Structureless_CenterMove(void);

double flatforce(double ); //checked
double Aduattractive(double , double );
double Adurepulsive(double , double );

void HemisphereInitContinue(void);
double d0=0.92;//solid lattice length
double zCut=0.3;
double Lgs6=powl(Lgs,6.0);
double Lgs12=powl(Lgs,12.0);

void timecounting(int);
double DCut=10.0;
double tol=0.1;

double AForce(double, double, double, double);
double AFUatt(double, double);
double AFUrep(double , double);
double FUatt(double, double);
double FUrep(double, double);





double ContactAngleCalculation(void);
void PrintContactAngle(double);
double contact_angle[100]; //to dynamically save contact angle 
int nCA; // the number of contact angle calculations 
Mol *fitmol,*Uppmol,*Fitatoms,*CentralDroplet,*LayerAtoms;

double circle[3], dcircle[3];
#define NFitPoints 100 //using 100 points to fit interface;
double InterfaceData[NFitPoints][2];
void circlefit(void);
void PrintInterface(void);

int aldle(void);

double ATA[3][3],ATB[3];//to solve ATA*x=ATB
int NCA=0;

double totalU;
double SteelePotential(double );
double AsymptoticPotential(double ,double ,double );
double FullPotential(double , double , double );
double GroovePotentialAtt(double ,double);
double GroovePotentialRep(double,double );
void PrintTotalEnergy();
void NoseHooverThermostat(void);
double gama=0;//
double Q=1;//coupling strength
void MomentumConservation(void);
void  PrintVelocityDistribution(void);
void VelocityRescaling(void);
void forceztotal(void);

double zforceArtificalStep(double ,double, double);
double yforceArtificalStep(double ,double, double);
int hstepGhost=int(DCut/(StepGap+StepWidth))+1; //=the half step ghost
void CountDensity(void);
void PrintDensity(void);
double DensityGap1;
double RoundDensity;
double RoundCenter;
double roundMaker(double );
void PrintGapCoord(int);
Mol *halfmolGap;
VecR SumCenter;
double ZCenterofMass;
void CenterOfMassCalc(void);
void PrintCenterofMass(void);
void velocity_Verlet(int );
void PrintFinish(void );
int  stepDesired2=21474836;
int  stepDesired3=21474836;
int  stepDesired4=21474836;
int  stepDesired5=21474836;
int  stepDesired6=21474836;
void PrintStructurelessCoordEvalution3(void );
void PrintStructurelessCoordEvalution4(void );
void PrintStructurelessCoordEvalution5(void );
void PrintStructurelessCoordEvalution6(void );
