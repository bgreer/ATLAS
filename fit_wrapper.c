#include "header.h"
#include "guess_table_16.h"

void free_output (outputs *outp);

/* need to pass:
 * guess file name
 * verbosity
 * input spectrum
 * output array
 */
// this function is called by the FORTRAN part of ATLAS
int fit_wrapper_ (double *spec_f, int *ntheta_f, int *nk_f, int *nnu_f, double *delta_k_f, double *delta_nu_f, 
		char *fname_guess, int *verbose_f, int *nmodes, float *dataout, int *narr, int *karr, int *kstart, int *kend)
{
	FILE *fp;
	int retval, ii, ij, junk, printback, verbose, ntheta, nk, nnu, ind;
	double delta_k, delta_nu, ***spec;
	int nummodes, *guess_k;
	double *guess_freq, *guess_amp, *guess_width;
	cl_tag tag;

	// file names
	char fname_pspec[FNAMELEN], fname_output[FNAMELEN];
	char fname_back[FNAMELEN];

	inputs inp; // struct containing all of the inputs for mrf
	outputs outp; // struct containing all of the outputs

	verbose = *verbose_f;
	ntheta = *ntheta_f;
	nk = *nk_f;
	nnu = *nnu_f;
	delta_k = *delta_k_f;
	delta_nu = *delta_nu_f;

	// make room for spectrum
	ind = 0;
	spec = (double***) malloc(nnu*sizeof(double**));
	for (ii=0; ii<nnu; ii++)
	{
		spec[ii] = (double**) malloc(nk*sizeof(double*));
		for (ij=0; ij<nk; ij++)
		{
			spec[ii][ij] = (double*) malloc(ntheta*sizeof(double));
			if (spec[ii][ij]==NULL)
			{
				printf("Error allocating memory.\n");
				return EXIT_FAILURE;
			} else {
				memcpy(spec[ii][ij], &(spec_f[ind]), ntheta*sizeof(double));
				ind += ntheta;
			}
		}
	}

	
	// load power spectrum data into input struct
	inp.pspec = spec;
	inp.pspec_ntheta = ntheta;
	inp.pspec_nk = nk;
	inp.pspec_nnu = nnu;
	inp.pspec_delta_k = delta_k;
	inp.pspec_delta_nu = delta_nu;
	
	// tell code that it needs to convert a pspec in log-power, kx-ky to one of linear-power, tht-k
	inp.convert_cart_polar = 0;

	// read a guess table
	//retval = read_model_file_testing ("16_deg.model", &nummodes, &guess_k, &guess_freq, &guess_amp, &guess_width, verbose);
  nummodes = guess_table_num;
  guess_k = (int*) malloc(nummodes * sizeof(int));
  guess_freq = (double*) malloc(nummodes * sizeof(double));
  guess_amp = (double*) malloc(nummodes * sizeof(double));
  guess_width = (double*) malloc(nummodes * sizeof(double));
  memcpy(guess_k, guess_table_k, nummodes * sizeof(int));
  memcpy(guess_freq, guess_table_freq, nummodes * sizeof(double));
  memcpy(guess_amp, guess_table_amp, nummodes * sizeof(double));
  memcpy(guess_width, guess_table_width, nummodes * sizeof(double));

	// load guess table into input struct
	inp.guess_nummodes = nummodes;
	inp.guess_k = guess_k; // wavenumber index
	inp.guess_freq = guess_freq;
	inp.guess_amp = guess_amp;
	inp.guess_width = guess_width;

	// set up run parameters
	// check meanings in interface.h
	inp.verbose = 0;//verbose;
	inp.dofits = 1;
	inp.fit_above = 2;
	inp.ac_cutoff = 4500.0;
	inp.detect_ridge = 1.;
	inp.kstart = *kstart;
	inp.kend = *kend;
	inp.ftol = 1e-8;
	inp.xtol = 1e-4;
	inp.gtol = 1e-8;
	inp.niter = 400;
	inp.output_covar = 0;

	// recording the entire covariance matrix is a pain
	// mostly from a data organization perspective
	// call fitting procedure
	if (verbose) printf("Calling fitting.. (%d-%d)\n", *kstart, *kend);
	mrf_fit(&inp, &outp);
	if (verbose) printf("Done fitting.\n");

	// do something with outputs
	if (verbose) printf("Total Modes Fit between (%d-%d) = %d\n", *kstart, *kend, outp.nummodes);

	*nmodes = outp.nummodes;

	for (ii=0; ii<outp.nummodes; ii++)
	{
		(narr)[ii] = outp.mode_n[ii];
		(karr)[ii] = outp.mode_k[ii];
		(dataout)[ii*4+0] = outp.mode_ux[ii];
		(dataout)[ii*4+1] = outp.mode_err_ux[ii];
		(dataout)[ii*4+2] = outp.mode_uy[ii];
		(dataout)[ii*4+3] = outp.mode_err_uy[ii];
	}

	/*fp = fopen(fname_output, "w");
	for (ii=0; ii<outp.nummodes; ii++)
	{
		// obviously there are more parameters that can be printed
		// check interface.h, outputs struct for all of the mode_* parameters
		fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%d\n", 
				outp.mode_k[ii], 
				(outp.mode_k[ii]+1)*delta_k, // the wavenumber in Mm^-1
				outp.mode_freq[ii] / ((outp.mode_k[ii]+1)*delta_k),
				outp.mode_freq[ii], outp.mode_err_freq[ii], 
				outp.mode_amp[ii], outp.mode_err_amp[ii], 
				outp.mode_width[ii], outp.mode_err_width[ii],
				outp.mode_ux[ii], outp.mode_err_ux[ii], 
				outp.mode_uy[ii], outp.mode_err_uy[ii], 
				outp.mode_fanis[ii], outp.mode_err_fanis[ii], 
				outp.mode_thtanis[ii], outp.mode_err_thtanis[ii], 
				0.0, 0.0, 
				outp.mode_n[ii],
				outp.mode_mpreturn[ii], // return code from optimization. used for data culling
				outp.mode_failure[ii]); // failure mode tag. also used for culling
	}
	fclose(fp);
*/

	// free memory? probably important..
	// TODO: free memory...

  free_output(&outp);

  free(guess_k);
  free(guess_freq);
  free(guess_amp);
  free(guess_width);
	return 0;
}

void free_output (outputs *outp)
{
  if (outp->mode_k) free(outp->mode_k);
  if (outp->mode_n) free(outp->mode_n);
  if (outp->mode_freq) free(outp->mode_freq);
  if (outp->mode_err_freq) free(outp->mode_err_freq);
  if (outp->mode_width) free(outp->mode_width);
  if (outp->mode_err_width) free(outp->mode_err_width);
  if (outp->mode_amp) free(outp->mode_amp);
  if (outp->mode_err_amp) free(outp->mode_err_amp);
  if (outp->mode_ux) free(outp->mode_ux);
  if (outp->mode_err_ux) free(outp->mode_err_ux);
  if (outp->mode_uy) free(outp->mode_uy);
  if (outp->mode_err_uy) free(outp->mode_err_uy);
  if (outp->mode_fanis) free(outp->mode_fanis);
  if (outp->mode_err_fanis) free(outp->mode_err_fanis);
  if (outp->mode_thtanis) free(outp->mode_thtanis);
  if (outp->mode_err_thtanis) free(outp->mode_err_thtanis);
  if (outp->mode_mpreturn) free(outp->mode_mpreturn);
  if (outp->mode_failure) free(outp->mode_failure);
  if (outp->back_k) free(outp->back_k);
  if (outp->back_amp) free(outp->back_amp);
  if (outp->back_err_amp) free(outp->back_err_amp);
  if (outp->back_cutoff) free(outp->back_cutoff);
  if (outp->back_err_cutoff) free(outp->back_err_cutoff);
  if (outp->back_index) free(outp->back_err_index);
  if (outp->back_fanis) free(outp->back_fanis);
  if (outp->back_err_fanis) free(outp->back_err_fanis);
  if (outp->back_thtanis) free(outp->back_thtanis);
  if (outp->back_err_thtanis) free(outp->back_err_thtanis);
}

int read_model_file_testing (char *fname, int *nummodes, int **guess_k, 
		double **guess_freq, double **guess_amp, double **guess_width, int verbose)
{
	FILE *fpmodel;
	int ii, ij, ik, il, num;
	float readfreq, readamp, readwidth, readk;

	if (verbose) printf("Reading model from %s\n", fname);
	num = 0;
	/* Run through once to count modes */
	fpmodel = fopen(fname, "r");
	if (fpmodel==NULL)
	{
		printf("ERROR: could not open model file: %s\n", fname);
		return EXIT_FAILURE;
	}
	while (fscanf(fpmodel, "%d\t%e\f%e\t%e\t%e\n", &ii, &readk, &readfreq, &readamp, &readwidth) != EOF)
		num++;
	rewind(fpmodel);

	/* allocate enough for all ridges, then load them */
	*nummodes = num;
	*guess_k = (int*) malloc(num*sizeof(int));
	*guess_freq = (double*) malloc(num*sizeof(double));
	*guess_amp = (double*) malloc(num*sizeof(double));
	*guess_width = (double*) malloc(num*sizeof(double));

	for (ii=0; ii<num; ii++)
	{
		il = fscanf(fpmodel, "%d\t%e\t%e\t%e\t%e\n", &ik, &readk, &readfreq, &readamp, &readwidth);
		if (readfreq<0.0 || readamp<0.0 || readwidth<0.0)
		{
			printf("ERROR: Invalid parameters in model, line=%d, k=%d\n", ii, ik);
			return EXIT_FAILURE;
		}
		(*guess_k)[ii] = ik;
		(*guess_freq)[ii] = readfreq;
		(*guess_amp)[ii] = readamp;
		(*guess_width)[ii] = readwidth;
	}
	fclose(fpmodel);
	return 0;
}

