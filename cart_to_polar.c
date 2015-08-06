#include <math.h>
#include <fftw3.h>
#include "header.h"

void convert_spectrum (inputs *inp, double ****spec, 
		int *ntheta, int *nk, int *nnu, double *delta_k, double *delta_nu)
{
	int ii, ij, ik, ntheta_full, ind;
	double xc, yc, x, y, dtheta, tht;
	double *line, *line2;
	fftw_complex *in, *out;
	fftw_plan p, p2;
	
	// copy dimensions that won't change
	*nnu = inp->pspec_nnu;
	*nk = inp->pspec_nk/2;
	*delta_k = inp->pspec_delta_k;
	*delta_nu = inp->pspec_delta_nu;

	// inp->pspec is allocated to [inp->pspec_nnu][inp->pspec_nk][inp->pspec_nk]
	// need to transform to polar coords
	
	// STEPS:
	// 1 - cart to polar
	// 2 - fourier subsample
	// 3 - fourier filter

	// pick number of points in theta to interpolate on to
	ntheta_full = 4*(*nk)/3;
	*ntheta = imax(ntheta_full/2, 16);
	if (inp->verbose) printf("Using %d points in theta\n", *ntheta);

	// allocate space
	*spec = (double***) malloc((*nnu) * sizeof(double**));
	for (ii=0; ii<*nnu; ii++)
	{
		(*spec)[ii] = (double**) malloc((*nk) * sizeof(double*));
		for (ij=0; ij<*nk; ij++)
		{
			(*spec)[ii][ij] = (double*) malloc((*ntheta) * sizeof(double));
		}
	}
	line = (double*) malloc(ntheta_full * sizeof(double));
	in = (fftw_complex*) malloc(ntheta_full * sizeof(fftw_complex));
	out = (fftw_complex*) malloc(ntheta_full * sizeof(fftw_complex));

	// center of spectrum
	xc = *nk;
	yc = *nk;
	dtheta = TWOPI / ntheta_full;

	// interpolate
	if (inp->verbose) printf("Unwrapping..");
	if (inp->verbose) fflush(stdout);
	p = fftw_plan_dft_1d(ntheta_full, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(ntheta_full, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (ii=0; ii<(*nnu); ii++)
	{
		for (ij=0; ij<(*nk); ij++)
		{
			// gather a line along theta
			for (ik=0; ik<ntheta_full; ik++)
			{
				tht = ik * dtheta;
				x = xc + ((double)ij+1.0)*cos(tht);
				y = yc + ((double)ij+1.0)*sin(tht);
				line[ik] = exp(interp(x, y, inp->pspec[ii], inp->pspec_nk));
			}
			// subsample that line
			// load line[] into in[]
			for (ik=0; ik<ntheta_full; ik++)
			{
				in[ik][0] = line[ik];
				in[ik][1] = 0.0;
			}
			// do fft
			fftw_execute(p);
			// mask
			// out[0] is 0th element
			for (ik=0; ik<(*ntheta - ntheta_full)*2; ik++)
			{
				ind = ik + ntheta_full/2 - (*ntheta - ntheta_full);
				out[ind][0] = 0.0;
				out[ind][1] = 0.0;
			}
			// do inverse fft
			fftw_execute(p2);
			// put result in (*spec)
			for (ik=0; ik<(*ntheta); ik++)
				(*spec)[ii][ij][ik] = in[ik*ntheta_full/(*ntheta)][0];
		}
	}
	fftw_destroy_plan(p);
	fftw_destroy_plan(p2);
		// debug
		 /*
		for (ij=0; ij<(*nnu); ij++)
		{
			for (ik=0; ik<(*ntheta); ik++)
			{
				printf("%d\t%d\t%f\n", ij, ik, (*spec)[ij][25][ik]);
			}
			printf("\n");
		}
		exit(-1);
*/

	if (inp->verbose) printf(" Done.\n");

	// free memory for old spectrum?
	free(line);
	free(in);
	free(out);
	
	// copy updated nk value back to input struct just in case
	inp->pspec_nk = *nk;
}

/* bicubic interpolation */
double interp (double x, double y, double **data, int dim)
{
	int ii, ij;
	int ix[4], iy[4];
	double r[4], weight, temp, res;

	if (x < 0.0 || x >= dim || y < 0.0 || y >= dim) return 0.0;

	// find nearest pixels
	ix[1] = imax(floor(x), 0);
	iy[1] = imax(floor(y), 0);
	ix[2] = imin(ix[1]+1, dim-1);
	iy[2] = imin(iy[1]+1, dim-1);
	ix[0] = imax(ix[1]-1,0);
	iy[0] = imax(iy[1]-1,0);
	ix[3] = imin(ix[1]+2,dim-1);
	iy[3] = imin(iy[1]+2,dim-1);


	// interpolate at current y
	for (ii=0; ii<4; ii++) // x index
	{
		r[ii] = 0.0;
		weight = 0.0;
		for (ij=0; ij<4; ij++) // y index
		{
			temp = kernel(iy[ij] - y);
			weight += temp;
			r[ii] += temp * data[ix[ii]][iy[ij]];
		}
		r[ii] /= weight; // normalize
	}

	// interpolate at current x
	for (ii=0; ii<4; ii++)
	{
		temp = kernel(ix[ii]-x);
		weight += temp;
		res += temp * r[ii];
	}

	return res / weight;
}

/* for bicubic interpolation */
double kernel (double x0)
{
	double x;
	x = fabs(x0);
	if (x <= 1.0)
		return 1.0 + x*x*(-2.5 + 1.5*x);
	else if (x < 2.0)
		return 2.0 + x*(-4.0 + x*(2.5 - 0.5*x));
	else
		return 0.0;
}

int imin (int a, int b)
{
	return (a > b) ? b : a;
}

int imax (int a, int b)
{
	return (a > b) ? a : b;
}
