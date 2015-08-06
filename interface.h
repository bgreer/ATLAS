
typedef struct
{
	// need to handle giving power spectrum, guess table, settings
	
	// power spectrum
	int pspec_ntheta, pspec_nk, pspec_nnu;
	double pspec_delta_k, pspec_delta_nu;
	double ***pspec; // pspec[nnu][nk][nk]

	// flag to say if the power spectrum needs to
	//  be converted from cartesian to polar coords
	int convert_cart_polar;

	// guess table
	int guess_nummodes; // total number of modes to fit
	int *guess_k; // wavenumber bin
	double *guess_freq, *guess_amp, *guess_width;

	// other settings
	
	// verbose mode
	int verbose;
	
	// perform optimization, 0=don't, 1=do
	// useful for testing guess table
	int dofits;
	
	// fit peaks above acoustic cutoff frequency
	// 0 = don't fit anything above
	// 1 = simple lorentzians (force velocity = 0)
	// 2 = full peaks everywhere
	int fit_above;
	
	// acoustic cutoff frequency, relevant for above setting
	// units = uHz
	double ac_cutoff;
	
	// Auto-detect ridge number using table
	// 0 = Don't, every n = 0
	// 1 = Do
	int detect_ridge;
	
	// which k-bin to start fitting at
	int kstart;
	
	// which k-bin to stop fitting at
	int kend;
	
	// optimization tolerances (likelihood, parameter, orthogonality)
	double ftol, xtol, gtol;
	
	// maximum number of iterations for optimization
	int niter;

	// output covariance matrices
	int output_covar;

} inputs;

typedef struct
{
	// need to handle returning fitted modes, background profiles, covariance matrices
	int nummodes;
	// next are 1d arrays of size >= nummodes
	// will be allocated by the fitting code
	int *mode_k;
	int *mode_n;
	double *mode_freq, *mode_err_freq;
	double *mode_amp, *mode_err_amp;
	double *mode_width, *mode_err_width;
	double *mode_ux, *mode_err_ux;
	double *mode_uy, *mode_err_uy;
	double *mode_fanis, *mode_err_fanis;
	double *mode_thtanis, *mode_err_thtanis;
	int *mode_mpreturn, *mode_failure;

	// background parameters
	int numback;
	// following are arrays of size >= numback
	int *back_k;
	double *back_amp, *back_err_amp;
	double *back_cutoff, *back_err_cutoff;
	double *back_index, *back_err_index;
	double *back_fanis, *back_err_fanis;
	double *back_thtanis, *back_err_thtanis;

	// covariance matrices
	
} outputs;
