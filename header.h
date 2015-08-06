#include "/shared/cfitsio-3.24/include/fitsio.h"
#include "mpfit.h"
#include "interface.h"

/* For later, optimizing */
/* NOT ACTUALLY USED. YET. */
/*
#define NNU (1152)
#define NK (192)
#define NTHETA (20)
*/

#define NPEAK (7) /* number of parameters for each peak */
#define NBACK (8) /* number of parameters for background terms */

/* Set order of parameters in parameter file */
#define PARAM_SPECTRUM (0)
#define PARAM_MODEL (1)
#define PARAM_OUTPUT (2)
#define PARAM_SILENT (3)
#define PARAM_KRANGE (4)
#define PARAM_FTOL (5)
#define PARAM_XTOL (6)
#define PARAM_GTOL (7)
#define PARAM_NITER (8)
#define PARAM_DEBUG (9)
#define PARAM_COVAR (10)
#define PARAM_BACK (11)
#define PARAM_FIT (12)
#define PARAM_CUTOFF (13)
#define PARAM_FITABOVE (14)
#define PARAM_RIDGE (15)

#define PI    (3.14159265358979324)
#define TWOPI (6.28318530717958648)

#define FNAMELEN 512

/* Types of tags that can be parsed */
#define TAGTYPE_BOOL 0
#define TAGTYPE_INT 1
#define TAGTYPE_FLOAT 2
#define TAGTYPE_STRING 3

/* Return values for error */
#define PARSE_ERROR -1
#define TAGERR_NAME -2
#define TAGERR_TYPE -3
#define TAGERR_DATA -4

/* used for reading command-line args */
typedef struct
{
	char *name; /* text to look for */
	int type; /* data type, see above #defines */
	void *data; /* the actual data */
} cl_tag;


struct kslice
{
	double **data;
	int start, end, ntheta, n;
	double delta_nu, k;
	struct params* par;
};

struct params
{
	int verbose, dofits, fit_above, detect_ridge;
	double ac_cutoff, multfac;
	int kstart, kend;
	double ftol, xtol, gtol;
	int niter;
	int output_covar;
};

void usage (char *argv);

/* mrf_fit.c */
int mrf_fit (inputs *inp, outputs *outp);
void mpreturn_translate (int mpreturn);
void normalize(double ***pol, struct params *par, int ntheta, int nk, int nnu);

/* cart_to_polar.c */
void convert_spectrum (inputs *inp, double ****spec, 
		int *ntheta, int *nk, int *nnu, double *delta_k, double *delta_nu);
double interp (double x, double y, double **data, int dim);
double kernel (double x0);
int imin (int a, int b);
int imax (int a, int b);

/* fit.c */
int fit_peak (struct params* p, double *freq, double *amp, double *width, 
	double*** pol, double delta_nu, double delta_k, int nnu, int ntheta, int k);
int funk_single(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int fit_back (struct params* p, double* amp, double* cutoff, double* power, 
	double*** pol, double delta_nu, int nnu, int ntheta, int k);
int funk_back(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* io.c */
void output_debug (struct params* p, double***pol, int ntheta, 
	int nk, int nnu, int k, int m, int n, double* x, double delta_nu, double delta_k);
void output_matrix (double* covar, int n, struct params* p, int ident);
int read_model_input (inputs *inp, struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width);
/*  depreciated: */
int read_model_file (struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width);
void read_fits_input (inputs *inp, struct params *p, double ****spec, 
		int *ntheta, int *nk, int *nnu, double *delta_k, double *delta_nu);
/*  depreciated: */
int read_fits_file (double**** spec, struct params* p, 
		int* ntheta, int* nk, int* nnu, double* delta_k, double* delta_nu);
void read_param_input (inputs *inp, struct params*p);
/*  depreciated: */
void read_param_file (char* fname, struct params* p);
void trim (char* str);

/* function.c */
int funk(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* errbars.c */
void errbars (int numparams, double *p, struct kslice *ks, double **covar);
void fisher (int numparams, double *p, struct kslice *ks, double *covar);
double model (double *p, struct kslice *ks, int inu, int itht);
void dm (double *p, struct kslice *ks, int inu, int itht, double *d);
void d2m (double *p, struct kslice *ks, int inu, int itht, double **d);

/* ridge.c */
int detect_ridge (double k, double nu);
double distance_to_ridge (double k, double nu, double a, double b, double c);


/* parse_input.c */
void set_parameters (int argc, char *argv[], struct params *par);
int parse_params (int argc, char **argv, int numtags, cl_tag* t);
int parse_params_action (char **argv, int ii, cl_tag t);
