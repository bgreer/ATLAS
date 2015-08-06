#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"


int mrf_fit (inputs *inp, outputs *outp)
{
	struct params par;
	int status = 0;
	int ii, ij, ik;
	int ntheta, nk, nnu;
	double delta_nu, delta_k;
	double ***pol;
	double **freq, **amp, **width;
	int *numridges, **fit_type, *ridgenum;

	int mpreturn, index;
	int failure_mode;
	struct kslice subsection;
	double* param;
	mp_result *mpres;
	mp_config *mpconf;
	mp_par *bounds;
	double *xerror;
	double *covar;
	mpres = malloc(sizeof(mp_result));
	mpconf = malloc(sizeof(mp_config));

	memset(mpres, 0x00, sizeof(mp_result));
	memset(mpconf, 0x00, sizeof(mp_config));

	int output_counter, output_back_counter;

	/* Copy settings data from input struct */
	read_param_input(inp, &par);


	/* Copy spectrum data from input struct */
	if (inp->convert_cart_polar)
		convert_spectrum (inp, &pol, &ntheta, &nk, &nnu, &delta_k, &delta_nu);
	else
		read_fits_input(inp, &par, &pol, &ntheta, &nk, &nnu, &delta_k, &delta_nu);

	/* normalize the power spectrum to avoid numerical issues */
	normalize(pol, &par, ntheta, nk, nnu);
	
	/* Copy guess table from input struct */
	read_model_input(inp, &par, nk, &numridges, &freq, &amp, &width);

	/* Prep the output struct */
	output_counter = 0;
	outp->nummodes = 0;
	outp->mode_k = (int*) malloc(inp->guess_nummodes * sizeof(int));
	outp->mode_n = (int*) malloc(inp->guess_nummodes * sizeof(int));
	outp->mode_freq = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_freq = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_amp = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_amp = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_width = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_width = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_ux = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_ux = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_uy = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_uy = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_fanis = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_fanis = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_thtanis = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_err_thtanis = (double*) malloc(inp->guess_nummodes * sizeof(double));
	outp->mode_mpreturn = (int*) malloc(inp->guess_nummodes * sizeof(int));
	outp->mode_failure = (int*) malloc(inp->guess_nummodes * sizeof(int));
	// background stuff
	output_back_counter = 0;
	outp->numback = 0;
	outp->back_k = (int*) malloc(inp->pspec_nk * sizeof(int));
	outp->back_amp = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_err_amp = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_cutoff = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_err_cutoff = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_index = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_err_index = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_fanis = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_err_fanis = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_thtanis = (double*) malloc(inp->pspec_nk * sizeof(double));
	outp->back_err_thtanis = (double*) malloc(inp->pspec_nk * sizeof(double));

	/* Check kstart and kend */
	if (par.kstart > par.kend) par.kend = par.kstart;
	if (par.kstart < 0) par.kstart = 0;
	if (par.kstart >= nk) par.kstart = nk-1;
	if (par.kstart < 0) par.kend = 0;
	if (par.kend >= nk) par.kend = nk-1;

	/* Decide which peaks to fit in what way */
	fit_type = malloc(nk*sizeof(int*));
	for (ii=0; ii<nk; ii++)
	{
		fit_type[ii] = 0;
		if (numridges[ii] > 0)
		{
			fit_type[ii] = malloc(numridges[ii]*sizeof(int));
			for (ij=0; ij<numridges[ii]; ij++)
			{
				if (par.fit_above==2)
					fit_type[ii][ij] = 1; /* Full line profile */
				else if (par.fit_above==1 && freq[ii][ij] >= par.ac_cutoff)
					fit_type[ii][ij] = 2; /* Simple line profile */
				else if (par.fit_above==0 && freq[ii][ij] >= par.ac_cutoff)
					fit_type[ii][ij] = 0; /* Do not fit */
				else
					fit_type[ii][ij] = 1;
			}
		}
	}

	if (par.verbose) printf("\nBeginning optimization.\n");


/* BEGIN LOOP OVER K-BINS */

	for (ij=par.kstart; ij<=par.kend; ij++)
	{
		// numridges[ij] gives the number of modes to fit at the current k-bin (ij)
		if (numridges[ij]>0)
		{
			if (par.verbose) printf("Starting fit of %d ridges at k=%d\n", numridges[ij], ij);
			param = (double*) malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			bounds = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(mp_par));
			xerror = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			ridgenum = malloc(numridges[ij]*sizeof(int));

			/* Do rough fit of single peaks */
			if (par.verbose) printf("\tDoing single ridge estimates.\n");
			for (ii=0; ii<numridges[ij]; ii++)
			{
				if (par.dofits)
					fit_peak(&par, &(freq[ij][ii]), &(amp[ij][ii]), &(width[ij][ii]), pol, 
						delta_nu, delta_k, nnu, ntheta, ij);
			}
			
			/* Do rough fit of background */
			if (par.verbose) printf("\tFitting background at low frequency.\n");
			fit_back(&par, &(param[numridges[ij]*NPEAK]), &(param[numridges[ij]*NPEAK+1]), 
				&(param[numridges[ij]*NPEAK+2]), pol, delta_nu, nnu, ntheta, ij);

			/* Set multifit contraints */
			/* if a parameter is not needed it is perfectly fine to set .fixed=1 to pin it to the value you give */
			for (ii=0; ii<numridges[ij]; ii++)
			{
				/* NPEAK is declared in header.h, gives number of parameters per peak (mode) */
				/* for the current mode, set default fit parameters */
				/* otherwise, set .limited[0] and .limited[1] to specify if a lower or upper constraint exists */
				/* then set .limits[0] and/or .limits[1] to specify the constraint values */
				/* bounds[i] is a struct containing info for the i-th parameter */
				/* parameters for peak 0 come first, then peak 1.... then after peaks comes the background params */
				for (ik=0; ik<NPEAK; ik++)
				{
					bounds[ii*NPEAK+ik].fixed = 0;
					bounds[ii*NPEAK+ik].side = 3; /* 3=analytical derivs, 0=numerical derivs */
					bounds[ii*NPEAK+ik].deriv_debug = 0;
					bounds[ii*NPEAK+ik].relstep = 0.0;
					bounds[ii*NPEAK+ik].step = 0.0;
					bounds[ii*NPEAK+ik].deriv_reltol = 0.0;
					bounds[ii*NPEAK+ik].deriv_abstol = 0.0;
				}

				/* 0 - set frequency */
				param[ii*NPEAK] = freq[ij][ii]; /* from single-ridge fitting */
				bounds[ii*NPEAK].limited[0] = bounds[ii*NPEAK].limited[1] = 1;
				if (ii==0)
					bounds[ii*NPEAK].limits[0] = 0.0;
				else 
					bounds[ii*NPEAK].limits[0] = 0.5*(freq[ij][ii] + freq[ij][ii-1]);
				if (ii==numridges[ij]-1)
					bounds[ii*NPEAK].limits[1] = nnu*delta_nu;
				else
					bounds[ii*NPEAK].limits[1] = 0.5*(freq[ij][ii] + freq[ij][ii+1]);

				/* 1 - set amplitude */
				param[ii*NPEAK+1] = amp[ij][ii]; /* from single-ridge fitting */
				bounds[ii*NPEAK+1].limited[0] = 1;
				bounds[ii*NPEAK+1].limited[1] = 0;
				bounds[ii*NPEAK+1].limits[0] = 0.0; /* dont go below 0 */

				/* 2 - set width */
				param[ii*NPEAK+2] = width[ij][ii]; /* from single-ridge fitting */
				bounds[ii*NPEAK+2].limited[0] = bounds[ii*NPEAK+2].limited[1] = 1;
				bounds[ii*NPEAK+2].limits[0] = width[ij][ii]/3.0;
				bounds[ii*NPEAK+2].limits[1] = width[ij][ii]*3.0;

				/* 3,4 - set velocities */
				if (param[ii*NPEAK] < par.ac_cutoff || par.fit_above != 1)
				{
					param[ii*NPEAK+3] = 0.0;
					bounds[ii*NPEAK+3].limited[0] = bounds[ii*NPEAK+3].limited[1] = 1;
					bounds[ii*NPEAK+3].limits[0] = -1000.0;
					bounds[ii*NPEAK+3].limits[1] = 1000.0;
					bounds[ii*NPEAK+3].step = 5.0;
					param[ii*NPEAK+4] = 0.0;
					bounds[ii*NPEAK+4].limited[0] = bounds[ii*NPEAK+4].limited[1] = 1;
					bounds[ii*NPEAK+4].limits[0] = -1000.0;
					bounds[ii*NPEAK+4].limits[1] = 1000.0;
					bounds[ii*NPEAK+4].step = 5.0;
				} else {
					param[ii*NPEAK+3] = 0.;
					param[ii*NPEAK+4] = 0.;
					bounds[ii*NPEAK+3].fixed = 1;
					bounds[ii*NPEAK+4].fixed = 1;
					bounds[ii*NPEAK+3].limited[0] = bounds[ii*NPEAK+3].limited[1] = 0;
					bounds[ii*NPEAK+4].limited[0] = bounds[ii*NPEAK+4].limited[1] = 0;
				}
				/* 5,6 - set theta variation params */
				param[ii*NPEAK+5] = 0.1;
				bounds[ii*NPEAK+5].limited[0] = 1;
				bounds[ii*NPEAK+5].limited[1] = 1;
				bounds[ii*NPEAK+5].limits[0] = -0.999;
				bounds[ii*NPEAK+5].limits[1] = 0.999;
				bounds[ii*NPEAK+5].fixed = 0; /* this will pin the anisotropy to the init value if set to 1 */
				/* because the fraction is allowed to go negative and the direction is allowed
				 * to go anywhere, there is a surface of constant likelihood between these two
				 * It's probably ok. Maybe. */
				param[ii*NPEAK+6] = PI/2.0;
				bounds[ii*NPEAK+6].limited[0] = bounds[ii*NPEAK+6].limited[1] = 0;

				for (ik=0; ik<7; ik++)
					if (bounds[ii*NPEAK+ik].limited[0] && param[ii*NPEAK+ik] < bounds[ii*NPEAK+ik].limits[0])
						printf("(inconsistent bounds/init value) error on param %d of peak %d\t%e %e\n", ik, ii, param[ii*NPEAK+ik], bounds[ii*NPEAK+ik].limits[0]);
			}

			/* set constraints for background */
			/* NBACK is set in header.h, gives number of parameters for the background function */
			for (ii=0; ii<NBACK; ii++)
			{
				bounds[numridges[ij]*NPEAK+ii].fixed = 0;
				bounds[numridges[ij]*NPEAK+ii].side = 3; /* 3=analytical deriv, 0=numerical deriv */
				bounds[numridges[ij]*NPEAK+ii].deriv_debug = 0;
				bounds[numridges[ij]*NPEAK+ii].relstep = 0.0;
				bounds[numridges[ij]*NPEAK+ii].step = 0.0;
				bounds[numridges[ij]*NPEAK+ii].deriv_reltol = 0.0;
				bounds[numridges[ij]*NPEAK+ii].deriv_abstol = 0.0;
			}

			/* Values and bounds for background power law */
			/* power law amplitude */
			bounds[numridges[ij]*NPEAK].limited[0] = 1;
			bounds[numridges[ij]*NPEAK].limited[1] = 0;
			bounds[numridges[ij]*NPEAK].limits[0] = 0.0;
			/* power law cutoff frequency */
			bounds[numridges[ij]*NPEAK+1].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+1].limits[0] = 1e-8;
			if (param[numridges[ij]*NPEAK+1] < bounds[numridges[ij]*NPEAK+1].limits[0])
				param[numridges[ij]*NPEAK+1] = bounds[numridges[ij]*NPEAK+1].limits[0];
			bounds[numridges[ij]*NPEAK+1].limited[1] = 0;
			/* power law index */
			bounds[numridges[ij]*NPEAK+2].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+2].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+2].limits[0] = 1e-5;
			if (param[numridges[ij]*NPEAK+2] < bounds[numridges[ij]*NPEAK+2].limits[0])
				param[numridges[ij]*NPEAK+2] = bounds[numridges[ij]*NPEAK+2].limits[0];
			/* anisotropy fraction */
			param[numridges[ij]*NPEAK+3] = 0.01;
			bounds[numridges[ij]*NPEAK+3].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+3].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+3].limits[0] = -0.9999; /* aaaah */
			bounds[numridges[ij]*NPEAK+3].limits[1] = 0.9999;
			bounds[numridges[ij]*NPEAK+3].fixed = 0; /* this will pin background anis to 0 if set to 1 */
			/* anisotropy direction */
			param[numridges[ij]*NPEAK+4] = PI/2.0;
			bounds[numridges[ij]*NPEAK+4].limited[0] = bounds[numridges[ij]*NPEAK+4].limited[1] = 0;

			/* Values and bounds for background lorentzian */
			param[numridges[ij]*NPEAK+5] = 0.01*param[numridges[ij]*NPEAK];
			bounds[numridges[ij]*NPEAK+5].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+5].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+5].limits[0] = 0.0;
			bounds[numridges[ij]*NPEAK+5].limits[1] = 100.;
			param[numridges[ij]*NPEAK+6] = 1500.;
			bounds[numridges[ij]*NPEAK+6].limited[0] = bounds[numridges[ij]*NPEAK+6].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+6].limits[0] = 400.;
			bounds[numridges[ij]*NPEAK+6].limits[1] = 3500.;/*nnu*delta_nu;*/
			param[numridges[ij]*NPEAK+7] = 1000.;
			bounds[numridges[ij]*NPEAK+7].limited[0] = bounds[numridges[ij]*NPEAK+7].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+7].limits[0] = 300.0;
			bounds[numridges[ij]*NPEAK+7].limits[1] = 5000.;

			/* Set background lorenztian to 0 */
			
			param[numridges[ij]*NPEAK+5] = 0.0;
			bounds[numridges[ij]*NPEAK+5].fixed = 1;
			param[numridges[ij]*NPEAK+6] = 1500.;
			bounds[numridges[ij]*NPEAK+6].fixed = 1;
			param[numridges[ij]*NPEAK+7] = 1000.;
			bounds[numridges[ij]*NPEAK+7].fixed = 1;

			/* DONE with setting parameter constraints */

			subsection.par = &par;
			subsection.n = numridges[ij];
			
			/* Load klice struct for passing to function */
			subsection.start = 400./delta_nu; /* start fitting at 400uHz */
			/* testing? */
			/*
			subsection.start = (param[0]-param[2])/delta_nu;
			if (subsection.start < 0) subsection.start = 0;
			*/

			/* stop fitting after the last peak */
			subsection.end = (param[(numridges[ij]-1)*NPEAK]+1.*param[(numridges[ij]-1)*NPEAK+2])/delta_nu;
			if (subsection.end > nnu-1) subsection.end = nnu-1;

			subsection.delta_nu = delta_nu;
			subsection.ntheta = ntheta;
			subsection.k = (ij+1)*delta_k;
			subsection.data = malloc((subsection.end-subsection.start+1)*sizeof(double*));
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				subsection.data[ii-subsection.start] = malloc(ntheta*sizeof(double));
				for (ik=0; ik<ntheta; ik++)
				{
					subsection.data[ii-subsection.start][ik] = pol[ii][ij][ik];
				}
			}

			/* Set optimization parameters */
			mpconf->ftol = par.ftol;
			mpconf->xtol = par.xtol;
			mpconf->gtol = par.gtol;
			mpconf->covtol = 1e-14;
			mpconf->maxiter = par.niter;
			mpconf->nofinitecheck = 1;
			mpreturn = 1;
			mpres->xerror = 0;
			mpres->covar = 0;
			mpres->resid = 0;


/* Perform optimization */
			if (par.dofits && par.verbose) printf("\tDoing optimization...\n");
		
			if (par.dofits) mpreturn = mpfit(&funk, 
					ntheta*(subsection.end-subsection.start+1), 
					numridges[ij]*NPEAK+NBACK, 
					param, 
					bounds, 
					mpconf, 
					&subsection, 
					mpres);
			/* Fitting is done, time to do post-processing and error checking */


			/* print params and bounds for consistency check */
			if (mpreturn == MP_ERR_INITBOUNDS)
			{
				for (ii=0; ii<NBACK+NPEAK*numridges[ij]; ii++)
				{
					if ((param[ii] < bounds[ii].limits[0] && bounds[ii].limited[0]) || (param[ii] > bounds[ii].limits[1] && bounds[ii].limited[1]))
						printf("%d %e %e %e\n", ii, bounds[ii].limits[0], param[ii], bounds[ii].limits[1]);
				}
			}

			/* Check return value of mpfit */
			if (par.verbose) mpreturn_translate(mpreturn);
		
			/* Output optimization details */
			if (par.verbose) 
			{
				printf("\tStarting LH = %e\n", 
					mpres->orignorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tFinal LH = %e\n", 
					mpres->bestnorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tNumber of iterations = %d\n", mpres->niter);
				printf("\tNumber of function evals = %d\n", mpres->nfev);
				printf("\tNumber of pegged parameters = %d\n", mpres->npegged);
			}


			/* Determine which parameters are pegged */
			if (mpres->npegged)
			{
				for (ii=0; ii<numridges[ij]; ii++)
				{
					for (ik=0; ik<NPEAK; ik++)
					{
						if (bounds[ii*NPEAK+ik].limited[0] && param[ii*NPEAK+ik] == bounds[ii*NPEAK+ik].limits[0])
							if (par.verbose)
								printf("\tParam %d of peak %d pegged low\n", ik, ii);
						if (bounds[ii*NPEAK+ik].limited[1] && param[ii*NPEAK+ik] == bounds[ii*NPEAK+ik].limits[1])
							if (par.verbose)
								printf("\tParam %d of peak %d pegged high\n", ik, ii);
					}
				}
				for (ik=0; ik<NBACK; ik++)
				{
					if (bounds[numridges[ij]*NPEAK+ik].limited[0] && 
							param[numridges[ij]*NPEAK+ik] == bounds[numridges[ij]*NPEAK+ik].limits[0])
						if (par.verbose)
							printf("\tParam %d of back pegged low\n", ik);
					if (bounds[numridges[ij]*NPEAK+ik].limited[1] && 
							param[numridges[ij]*NPEAK+ik] == bounds[numridges[ij]*NPEAK+ik].limits[1])
						if (par.verbose)
							printf("\tParam %d of back pegged high\n", ik);
				}
			}

			param[numridges[ij]*NPEAK+5] += 1e-10;

			/* Compute error bars */
			errbars(numridges[ij]*NPEAK+NBACK, param, &subsection, &covar);
			
			/* Load covar into xerror because I'm lazy */
			for (ii=0; ii<numridges[ij]*NPEAK+NBACK; ii++)
				xerror[ii] = sqrt(fabs(covar[ii*(numridges[ij]*NPEAK+NBACK)+ii]));

			/* Post-process anisotropy parameters (negative amplitude -> pi/2 shift) */
			for (ii=0; ii<numridges[ij]; ii++)
			{
				if (param[ii*NPEAK+5] < 0.0)
				{
					param[ii*NPEAK+6] += PI/2.0;
					param[ii*NPEAK+5] = fabs(param[ii*NPEAK+5]);
				}
				param[ii*NPEAK+6] = fmod(param[ii*NPEAK+6], PI);
				if (param[ii*NPEAK+6] < 0.0)
					param[ii*NPEAK+6] += PI;
			}
			if (param[numridges[ij]*NPEAK+3] < 0.0)
			{
				param[numridges[ij]*NPEAK+4] += PI/2.0;
				param[numridges[ij]*NPEAK+3] = fabs(param[numridges[ij]*NPEAK+3]);
			}
			param[numridges[ij]*NPEAK+4] = fmod(param[numridges[ij]*NPEAK+4], PI);
			if (param[numridges[ij]*NPEAK+4] < 0.0)
				param[numridges[ij]*NPEAK+4] += PI;

			/* Auto-detect ridge */
			if (par.detect_ridge)
			{
				for (ii=0; ii<numridges[ij]; ii++)
					ridgenum[ii] = detect_ridge((ij+1)*delta_k, param[ii*NPEAK]);
			} else {
				for (ii=0; ii<numridges[ij]; ii++)
					ridgenum[ii] = 0;
			}

			/* Print fit debug */
			/* Not relevant in DRMS version */
//			output_debug(&par, pol, ntheta, nk, nnu, ij, ntheta*(subsection.end-subsection.start+1), 
//					numridges[ij]*NPEAK+NBACK, 
//					param, delta_nu, delta_k);
			

			/* Print covariance matrix */
			//if (par.covarfname) output_matrix(covar, numridges[ij]*NPEAK+NBACK, &par, ij);
			if (inp->output_covar)
			{
				// too involved for now
			}


			/* Output fit to file if valid */
			for (ii=0; ii<numridges[ij]; ii++)
			{
				failure_mode = 0;
				if (param[ii*NPEAK+2] < 7.5)
				{
					printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
					printf("Width too small (%f)\n", param[ii*NPEAK+2]);
					failure_mode = 1;
				} else if (param[ii*NPEAK+2] > 900.) {
					printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
					printf("Width too large (%f)\n", param[ii*NPEAK+2]);
					failure_mode = 2;
				} else if (abs(param[ii*NPEAK]-freq[ij][ii])/freq[ij][ii] >= 0.2) {
					printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
					printf("Frequency drifted too much from model\n");
					failure_mode = 3;
				}/* else if (xerror[ii*NPEAK]<=0.0 || xerror[ii*NPEAK+1]<=0.0 || 
						xerror[ii*NPEAK+2]<=0.0 || xerror[ii*NPEAK+3]<=0.0 || 
						xerror[ii*NPEAK+4]<=0.0) {
					printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
					printf("Zero error\n");
				}*/ else {
					/* place results in output struct */
					outp->mode_k[output_counter] = ij;
					outp->mode_n[output_counter] = ridgenum[ii];
						
					outp->mode_freq[output_counter] = param[ii*NPEAK];
					outp->mode_err_freq[output_counter] = xerror[ii*NPEAK];
						
					outp->mode_amp[output_counter] = param[ii*NPEAK+1];
					outp->mode_err_amp[output_counter] = xerror[ii*NPEAK+1];

					outp->mode_width[output_counter] = param[ii*NPEAK+2];
					outp->mode_err_width[output_counter] = xerror[ii*NPEAK+2];

					outp->mode_ux[output_counter] = param[ii*NPEAK+3];
					outp->mode_err_ux[output_counter] = xerror[ii*NPEAK+3];

					outp->mode_uy[output_counter] = param[ii*NPEAK+4];
					outp->mode_err_uy[output_counter] = xerror[ii*NPEAK+4];

					outp->mode_fanis[output_counter] = param[ii*NPEAK+5];
					outp->mode_err_fanis[output_counter] = xerror[ii*NPEAK+5];

					outp->mode_thtanis[output_counter] = param[ii*NPEAK+6];
					outp->mode_err_thtanis[output_counter] = xerror[ii*NPEAK+6];

					outp->mode_mpreturn[output_counter] = mpreturn;
					outp->mode_failure[output_counter] = failure_mode;
					output_counter++;
				}
			}
			

			/* Output background params to separate file */
			outp->back_k[output_back_counter] = ij;
			outp->back_amp[output_back_counter] = param[numridges[ij]*NPEAK+0];
			outp->back_err_amp[output_back_counter] = xerror[numridges[ij]*NPEAK+0];
			outp->back_cutoff[output_back_counter] = param[numridges[ij]*NPEAK+1];
			outp->back_err_cutoff[output_back_counter] = xerror[numridges[ij]*NPEAK+1];
			outp->back_index[output_back_counter] = param[numridges[ij]*NPEAK+2];
			outp->back_err_index[output_back_counter] = xerror[numridges[ij]*NPEAK+2];
			outp->back_fanis[output_back_counter] = param[numridges[ij]*NPEAK+3];
			outp->back_err_fanis[output_back_counter] = xerror[numridges[ij]*NPEAK+3];
			outp->back_thtanis[output_back_counter] = param[numridges[ij]*NPEAK+4];
			outp->back_err_thtanis[output_back_counter] = xerror[numridges[ij]*NPEAK+4];
			output_back_counter++;


			/* Free some memory for next k */
			free(param);
			free(bounds);
			free(xerror);
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				free(subsection.data[ii-subsection.start]);
			}
			free(subsection.data);
			free(covar);
			free(ridgenum);
		}
	}

/* END LOOP OVER K-BINS */

	outp->nummodes = output_counter;
	outp->numback = output_back_counter;
	
	/* Free ALL THE THINGS */
	for (ii=0; ii<nk; ii++)
	{
		if (numridges[ii]>0)
		{
			free(freq[ii]);
			free(amp[ii]);
			free(width[ii]);
			free(fit_type[ii]);
		}
	}
	free(fit_type);
	free(freq);
	free(amp);
	free(width);
	free(numridges);
	free(mpres);
	free(mpconf);
	for (ii=0; ii<nnu; ii++)
	{
		for (ij=0; ij<nk; ij++)
		{
			free(pol[ii][ij]);
		}
		free(pol[ii]);
	}
	free(pol);
	return EXIT_SUCCESS;
}

void mpreturn_translate (int mpreturn)
{
	switch (mpreturn)
	{
		case MP_OK_CHI:
			printf("\tMPFIT convergence in chi-square\n");break;
		case MP_OK_PAR:
			printf("\tMPFIT convergence in parameter value\n");break;
		case MP_OK_DIR:
			printf("\tMPFIT convergence in orthogonality\n");break;
		case MP_MAXITER:
			printf("\tMPFIT max iterations reached\n");break;
		case MP_FTOL:
			printf("\tMPFIT ftol criteria reached\n");break;
		case MP_XTOL:
			printf("\tMPFIT xtol criteria reached\n");break;
		case MP_GTOL:
			printf("\tMPFIT gtol criteria reached\n");break;

		case MP_ERR_INPUT:
			printf("\tMPFIT ERROR: input parameter error\n");break;
		case MP_ERR_NAN:
			printf("\tMPFIT ERROR: function returned nan\n");break;
		case MP_ERR_MEMORY:
			printf("\tMPFIT ERROR: memory allocation error\n");break;
		case MP_ERR_INITBOUNDS:
			printf("\tMPFIT ERROR: guesses inconsistent with bounds\n");break;
		case MP_ERR_BOUNDS:
			printf("\tMPFIT ERROR: inconsistent bounds\n");break;
		case MP_ERR_PARAM:
			printf("\tMPFIT ERROR: input parameter error\n");break;
		case MP_ERR_DOF:
			printf("\tMPFIT ERROR: too few degrees of freedom\n");break;

		default:
			printf("\tMPFIT return = %d\n", mpreturn);break;
	}
}

/* determine appropriate normalization factor and apply it */
void normalize(double ***pol, struct params *par, int ntheta, int nk, int nnu)
{
	int ii, ij, ik;
	double sum = 0.0;
	for (ii=0; ii<nnu; ii++)
	{
		for (ij=0; ij<nk; ij++)
		{
			for (ik=0; ik<ntheta; ik++)
			{
				sum += pol[ii][ij][ik];
			}
		}
	}
	sum /= (nnu*nk*ntheta);
	par->multfac = sum;
	if (par->verbose) printf("Normalizing to %e\n", sum);
	for (ii=0; ii<nnu; ii++)
	{
		for (ij=0; ij<nk; ij++)
		{
			for (ik=0; ik<ntheta; ik++)
			{
				pol[ii][ij][ik] /= sum;
			}
		}
	}
}
