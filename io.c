#include <string.h>
#include <math.h>
#include "header.h"

/* no longer used in drms version */

void output_debug (struct params* p, double ***pol, int ntheta, int nk, int nnu, int k, int m, int n, double* x, double delta_nu, double delta_k)
{
	FILE* fp;
	int ii, ik, ij, num;
	double back1, back2, fit, denom, sp, no;
	
	fp = fopen("debug", "w");
	for (ii=0; ii<nnu; ii++)
	{
		back1 = back2 = fit = sp = no = 0.0;
		for (ik=0; ik<ntheta; ik++)
		{
			back1 += x[n-8]*(1.0+x[n-5]*cos(2.0*(TWOPI*ik/ntheta - x[n-4])))/(1.0+pow(ii*delta_nu/(x[n-7]),x[n-6]));
			back2 += 0.5*x[n-3]*x[n-1]/((ii*delta_nu-x[n-2])*(ii*delta_nu-x[n-2]) + 0.25*x[n-1]*x[n-1]);
			for (ij=0; ij<(n-NBACK)/NPEAK; ij++)
			{
				denom = ii*delta_nu
					+ k*delta_k*(x[ij*NPEAK+3]*cos(ik*TWOPI/ntheta)
						+x[ij*NPEAK+4]*sin(ik*TWOPI/ntheta))/TWOPI
					- x[ij*NPEAK];
				fit += 0.5*x[ij*NPEAK+2]*x[ij*NPEAK+1] / 
					(denom*denom + 0.25*x[ij*NPEAK+2]*x[ij*NPEAK+2]);
			}
			sp += pol[ii][k][ik];
		}

		fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ik, ii*delta_nu, 
			sp, 
			0.0, 
			back1+back2+fit, 
			back1, 
			back2, 
			((back1+back2+fit)-sp)/(back1+back2+fit));
	}
	fclose(fp);
}



/* Prints square matrix to output file (Useful for covariance matrix) */
/* depreciated in drms version */
/*
void output_matrix (double* covar, int n, struct params* p, int ident)
{
	FILE *fp;
	int ii, ik;
	
	fp = fopen(p->covarfname, "a");
	if (fp==NULL)
	{
		printf("ERROR: Could not open covar output file %s\n", p->covarfname);
		return;
	}
	fprintf(fp, "%d\n", ident);
	for (ii=0; ii<n; ii++)
	{
		for (ik=0; ik<n; ik++)
		{
			fprintf(fp, "%d\t%d\t%e\n", ii, ik, covar[ii*n+ik]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
*/

void read_fits_input (inputs *inp, struct params *p, double ****spec, 
		int *ntheta, int *nk, int *nnu, double *delta_k, double *delta_nu)
{
	*spec = inp->pspec;
	*ntheta = inp->pspec_ntheta;
	*nnu = inp->pspec_nnu;
	*nk = inp->pspec_nk;
	*delta_k = inp->pspec_delta_k;
	*delta_nu = inp->pspec_delta_nu;
}


/* Read FITS spectrum, load header keys into variables, load data cube into spec */
/* depreciated in DRMS version. use read_fits_input() */
/*
int read_fits_file (double ****spec, struct params* p, 
		int* ntheta, int* nk, int* nnu, double* delta_k, double* delta_nu)
{
	fitsfile *fptr;
	int ii, ij, ik;
	int status = 0;
	long coords[3];
	float *buff, read;
	double sum;

	if (p->verbose) printf("Reading FITS file: %s\n", p->fitsfname);
	
	fits_open_file(&fptr, p->fitsfname, READONLY, &status);
	
	if (status)
	{
		printf("Error opening FITS file: %s\n", p->fitsfname);
		fits_report_error(stdout, status);
		return EXIT_FAILURE;
	}
	
	fits_read_key(fptr, TINT, "NAXIS1", ntheta, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", nk, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS3", nnu, NULL, &status);
	
	if (status)
	{
		printf("Could not find FITS dimensions.\n");
		printf("Make sure the following keys are set: NAXIS1, NAXIS2, NAXIS3\n");
		return EXIT_FAILURE;
	}
	
	if (p->verbose)
		printf("FITS dimensions:\n\ttheta:\t%d\n\tk:\t%d\n\tnu:\t%d\n", *ntheta, *nk, *nnu);
	
	fits_read_key(fptr, TFLOAT, "DELTA_K", &read, NULL, &status);
	*delta_k = read;
	fits_read_key(fptr, TFLOAT, "DELTA_NU", &read, NULL, &status);
	*delta_nu = read;
	
	if (status)
	{
		printf("Could not find FITS keys DELTA_NU and/or DELTA_K.\n");
		return EXIT_FAILURE;
	}
	if (p->verbose) printf("\tDELTA_NU:\t%f\n\tDELTA_K:\t%f\n", *delta_nu, *delta_k);

// Allocate memory for entire FITS data cube 
	if (p->verbose)
		printf("\nAllocating %ld bytes for data cube.\n", sizeof(double)*(*ntheta)*(*nk)*(*nnu));
	if (p->verbose)
		printf("\t With %ld bytes of overhead.\n", sizeof(double*)*(*nk)*(*nnu)+sizeof(double**)*(*nnu));

	*spec = (double***) malloc((*nnu)*sizeof(double**));
	for (ii=0; ii<(*nnu); ii++)
	{
		(*spec)[ii] = (double**) malloc((*nk)*sizeof(double*));
		for (ij=0; ij<(*nk); ij++)
		{
			(*spec)[ii][ij] = (double*) malloc((*ntheta)*sizeof(double));
			if ((*spec)[ii][ij]==NULL)
			{
				printf("Error allocating memory.\n");
				return EXIT_FAILURE;
			}
		}
	}
	buff = (float*) malloc((*ntheta)*sizeof(float));
// Read FITS file into array 
	coords[0] = 1L;
	if (p->verbose) printf("Reading data cube into memory.\n");
	for (coords[2]=1; coords[2]<=(*nnu); coords[2]++)
	{
		for (coords[1]=1; coords[1]<=(*nk); coords[1]++)
		{
			fits_read_pix(fptr, TFLOAT, coords, (*ntheta), NULL, buff, NULL, &status);
			for (ik=0; ik<(*ntheta); ik++)
			{
				if (buff[ik] < 0.0 || isnan(buff[ik])) buff[ik] = 0.0;
				(*spec)[coords[2]-1][coords[1]-1][ik] = buff[ik];
			}
		}
	}
	fits_close_file(fptr, &status);
	free(buff);
	return EXIT_SUCCESS;
}
*/

int read_model_input (inputs *inp, struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width)
{
	int ii, ij, ik, il;
	double readfreq, readamp, readwidth, readk;
	double ampsum;

	if (par->verbose) printf("Parsing guess table..\n");
	*numridges = (int*) malloc(nk*sizeof(int));
	*freq = (double**) malloc(nk*sizeof(double*));
	*amp = (double**) malloc(nk*sizeof(double*));
	*width = (double**) malloc(nk*sizeof(double*));
	for (ii=0; ii<nk; ii++)
		(*numridges)[ii] = 0;
	/* Run through once to count ridges */
	for (ij=0; ij<inp->guess_nummodes; ij++)
	{
		ii = inp->guess_k[ij];
		if (ii>=0 && ii<nk)
		{
			(*numridges)[ii]++;
		} else {
			printf("ERROR: k-value out of range in model: %d\n", ii);
			return EXIT_FAILURE;
		}
	}
	/* allocate enough for all ridges, then load them */
	il = 0;
	for (ii=0; ii<nk; ii++)
	{
		if ((*numridges)[ii]>0)
		{
			(*freq)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*amp)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*width)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
		}
		ampsum = 0.0;
		for (ij=0; ij<(*numridges)[ii]; ij++)
		{
			ik = inp->guess_k[il];
			readfreq = inp->guess_freq[il];
			readamp = inp->guess_amp[il];
			readwidth = inp->guess_width[il];

			if (readfreq<0.0 || readamp<0.0 || readwidth<0.0)
			{
				printf("ERROR: Invalid parameters in model, k=%d\n", ii);
				return EXIT_FAILURE;
			}
			(*freq)[ii][ij] = readfreq;
			(*amp)[ii][ij] = readamp;
			(*width)[ii][ij] = readwidth;
			ampsum += (*amp)[ii][ij];
			il ++;
		}
		for (ij=0; ij<(*numridges)[ii]; ij++)
			(*amp)[ii][ij] /= ampsum;
	}
	return 0;
}

/* Read model file */
/* depreciated in DRMS version */
/*
int read_model_file (struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width)
{
	FILE *fpmodel;
	int ii, ij, ik, il;
	float readfreq, readamp, readwidth, readk;

	if (par->verbose) printf("Reading model from %s\n", par->modelfname);
	*numridges = (int*) malloc(nk*sizeof(int));
	*freq = (double**) malloc(nk*sizeof(double*));
	*amp = (double**) malloc(nk*sizeof(double*));
	*width = (double**) malloc(nk*sizeof(double*));
	for (ii=0; ii<nk; ii++)
		(*numridges)[ii] = 0;
	// Run through once to count ridges 
	fpmodel = fopen(par->modelfname, "r");
	if (fpmodel==NULL)
	{
		printf("ERROR: could not open model file: %s\n", par->modelfname);
		return EXIT_FAILURE;
	}
	while (fscanf(fpmodel, "%d\t%e\f%e\t%e\t%e\n", &ii, &readk, &readfreq, &readamp, &readwidth) != EOF)
	{
		if (ii>=0 && ii<nk)
		{
			(*numridges)[ii]++;
		} else {
			printf("ERROR: k-value out of range in model: %d\n", ii);
			return EXIT_FAILURE;
		}
	}
	rewind(fpmodel);
	// allocate enough for all ridges, then load them 
	for (ii=0; ii<nk; ii++)
	{
		if ((*numridges)[ii]>0)
		{
			(*freq)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*amp)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*width)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
		}
		for (ij=0; ij<(*numridges)[ii]; ij++)
		{
			il = fscanf(fpmodel, "%d\t%e\t%e\t%e\t%e\n", &ik, &readk, &readfreq, &readamp, &readwidth);
			if (readfreq<0.0 || readamp<0.0 || readwidth<0.0)
			{
				printf("ERROR: Invalid parameters in model, k=%d\n", ii);
				return EXIT_FAILURE;
			}
			(*freq)[ii][ij] = readfreq;
			(*amp)[ii][ij] = readamp/par->multfac;
			(*width)[ii][ij] = readwidth;
		}
	}
	fclose(fpmodel);
	return 0;
}
*/


/* take an input struct and put the input parameter stuff
 * into the params struct. seems silly, but I'd rather do
 * this than recode the whole thing to use the input struct
 * instead of the param struct. */
void read_param_input (inputs *inp, struct params *p)
{
	p->verbose = inp->verbose;
	p->dofits = inp->dofits;
	p->fit_above = inp->fit_above;
	p->detect_ridge = inp->detect_ridge;
	p->ac_cutoff = inp->ac_cutoff;
	p->kstart = inp->kstart;
	p->kend = inp->kend;
	p->ftol = inp->ftol;
	p->xtol = inp->xtol;
	p->gtol = inp->gtol;
	p->niter = inp->niter;
	p->output_covar = inp->output_covar;
}



/* Reads parameter file fname and loads data into parameter struct *p */
/* this is depreciated for the DRMS version
 * use read_param_input instead.
 * included for sentimental purposes only*/
/*
void read_param_file (char* fname, struct params* p)
{
	FILE *fp;
	char buffer[200];
	int index, ii;

	fp = fopen(fname, "r");
	if (fp==NULL)
	{
		printf("Could not open parameter file: %s\n", fname);
		return;
	}

	// Allocate space for filenames 
	p->fitsfname = malloc(200*sizeof(char));
	p->modelfname = malloc(200*sizeof(char));
	p->outfname = malloc(200*sizeof(char));
	p->debugfname = 0;
	p->covarfname = 0;
	p->multfac = 1.0;

	// Scan through for non-comment and non-blank lines 
	index = 0;
	while (fgets(buffer, sizeof(buffer), fp) != NULL)
	{
		if (buffer[0] != '#' && buffer[0] != '\n')
		{
			trim(buffer);
			switch (index)
			{
				case PARAM_SPECTRUM:
					strcpy(p->fitsfname, buffer);
					break;
				case PARAM_MODEL:
					strcpy(p->modelfname, buffer);
					break;
				case PARAM_OUTPUT:
					strcpy(p->outfname, buffer);
					break;
				case PARAM_SILENT:
					p->verbose = atoi(buffer);
					break;
				case PARAM_KRANGE:
					// split string by putting null-terminator where the space is 
					ii=0;
					while (ii<strlen(buffer) && buffer[ii]!=' ')
						ii++;
					buffer[ii] = 0;
					p->kstart = atoi(buffer);
					p->kend = atoi(buffer+ii+1); // pointer arithmetic 
					break;
				case PARAM_FTOL:
					p->ftol = atof(buffer);
					break;
				case PARAM_XTOL:
					p->xtol = atof(buffer);
					break;
				case PARAM_GTOL:
					p->gtol = atof(buffer);
					break;
				case PARAM_NITER:
					p->niter = atoi(buffer);
					break;
				case PARAM_DEBUG:
					if (strcmp(buffer, "0")!=0)
					{
						p->debugfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->debugfname, buffer);
					}
					break;
				case PARAM_COVAR:
					if (strcmp(buffer, "0")!=0)
					{
						p->covarfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->covarfname, buffer);
					}
					break;
				case PARAM_BACK:
					if (strcmp(buffer, "0")!=0)
					{
						p->backfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->backfname, buffer);
					}
					break;
				case PARAM_FIT:
					p->dofits = atoi(buffer);
					break;
				case PARAM_CUTOFF:
					p->ac_cutoff = (double) atof(buffer);
					break;
				case PARAM_FITABOVE:
					p->fit_above = atoi(buffer);
					break;
				case PARAM_RIDGE:
					p->detect_ridge = atoi(buffer);
					break;
			}
			index++;
		}
	}
	fclose(fp);

	// If verbose, print run params (for sanity checks) 
	if (p->verbose)
	{
		printf("\nRun parameters:\n");
		printf("\tAcoustic Cutoff Frequency: %6.1f uHz\n", p->ac_cutoff);
		switch (p->fit_above)
		{
			case 0:
				printf("\tNo peaks fit above this.\n"); break;
			case 1:
				printf("\tSimple peaks fit above this.\n"); break;
			case 2:
				printf("\tFull peaks fit above this.\n"); break;
		}
		if (p->debugfname)
			printf("\tDebug filename base: %s\n", p->debugfname);
		if (p->covarfname)
			printf("\tCovar filename base: %s\n", p->covarfname);
		printf("\n");
	}
}
*/


