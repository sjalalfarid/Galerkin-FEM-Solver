#ifndef MATRIXOPERATIONS


//#include "Nrutil.h"
//#include "GeometryManipulator.h"

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
int** NewMatrixI(int);
int** NewMatrixI(int n, int m);
void DeleteMatrix(int, int**);
void DeleteMatrix(int,int,int**);
double** NewMatrix(int);
double** NewMatrix(int n, int m);
void DeleteMatrix(int, double**);
void DeleteMatrix(int,int,double**);
int SolveSymmetric(int , double** , double* );
int SymmetricInverse (int n, double **A, double **Ainv);
void Invert(double **orig, double **ans, int n);
void SolveAxB(double **inputLHS, double *inputRHS, int n);
void MatrixTimes(double **ans, double **m1, double **m2, int row1, int col1, int row2, int col2);
void MatrixTimes(double *ans, double *m1, double **m2, int row1, int col1, int row2, int col2);
void MatrixTimes(double *ans, double **m1, double *m2, int row1, int col1, int row2, int col2);
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
double pythag(double a, double b);
float pythag(float a, float b);
double Determinant4x4(double **m);
double det3x3(double  a1,double  a2,double  a3,double  b1,double  b2,double 
				b3,double  c1,double  c2,double  c3 );
double det3x3(double **a);
double det2x2(double a,double  b,double  c,double  d);
//double ind(GeoInfoStruct A, GeoInfoStruct B, GeoInfoStruct C, GeoInfoStruct D);
void tetrahedron_barycentric_3d ( double x1, double y1, double z1, double x2, 
  double y2, double z2, double x3, double y3, double z3, double x4, double y4, 
  double z4, double px, double py, double pz,
	double &r, double &s, double &t, double &l);
int dmat_solve ( int n, int nrhs, double a[] );





#define MATRIXOPERATIONS
#endif

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) 
static double dsqrarg; 
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg) 
static double dmaxarg1,dmaxarg2; 
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2)) 
static double dminarg1,dminarg2; 
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
(dminarg1) : (dminarg2)) 
static float maxarg1,maxarg2; 
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2)) 
static float minarg1,minarg2; 
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
(minarg1) : (minarg2)) 
static long lmaxarg1,lmaxarg2; 
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
(lmaxarg1) : (lmaxarg2)) 
static long lminarg1,lminarg2; 
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
(lminarg1) : (lminarg2)) 
static int imaxarg1,imaxarg2; 
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
(imaxarg1) : (imaxarg2)) 
static int iminarg1,iminarg2; 
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2)) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 
//#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */ 
void nrerror(char error_text[]); 
float *vector(long nl, long nh); 
int *ivector(long nl, long nh); 
unsigned char *cvector(long nl, long nh); 
unsigned long *lvector(long nl, long nh); 
double *dvector(long nl, long nh); 
float **matrix(long nrl, long nrh, long ncl, long nch); 
double **dmatrix(long nrl, long nrh, long ncl, long nch); 
int **imatrix(long nrl, long nrh, long ncl, long nch); 
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, 
long newrl, long newcl); 
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch); 
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh); 
void free_vector(float *v, long nl, long nh); 
void free_ivector(int *v, long nl, long nh); 
void free_cvector(unsigned char *v, long nl, long nh); 
void free_lvector(unsigned long *v, long nl, long nh); 
void free_dvector(double *v, long nl, long nh); 
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch); 
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch); 
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch); 
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch); 
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch); 
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, 
long ndl, long ndh); 
//#else /* ANSI */ 
/* traditional - K&R */ 
//void nrerror(); 
//float *vector(); 
 
//#endif /* ANSI */ 
#endif /* _NR_UTILS_H_ */ 

