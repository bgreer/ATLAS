#include "header.h"
#include "math.h"

/* Given a mode wavenumber and frequency, determing which ridge it is (f-mode, p1, p2, etc) */
int detect_ridge (double k, double nu)
{
	int ii, closest;
	double distance, temp;
	/* coefficients of polynomials in log-log space for each ridge */
	int known = 14;
	double c[] = {7.87, 0.48, 0.00, /* f-mode */
				  8.08, 0.48, 0.05, /* p1 */
				  8.26, 0.47, 0.05, /* p2 */
				  8.41, 0.45, 0.04, 
				  8.55, 0.45, 0.03, 
				  8.67, 0.47, 0.04, 
				  8.77, 0.48, 0.04, 
				  8.85, 0.48, 0.04, 
				  8.94, 0.49, 0.04, 
				  9.02, 0.50, 0.04, 
				  9.08, 0.50, 0.04,  /* p10 */
				  9.14, 0.51, 0.04,  /* p11 */
				  9.20, 0.51, 0.04,  /* p12 */
				  9.25, 0.51, 0.04   /* p13 */
				  };

	distance = 1e9;
	for (ii=0; ii<known; ii++)
	{
		if ((temp=distance_to_ridge(k, nu, c[ii*3], c[ii*3+1], c[ii*3+2])) < distance)
		{
			distance = temp;
			closest = ii;
		}
	}
	return closest;
}

double distance_to_ridge (double k, double nu, double a, double b, double c)
{
	double poly, logk;

	logk = log(k);
	poly = a + b*logk + c*logk*logk;

	return fabs(poly-log(nu));
}
