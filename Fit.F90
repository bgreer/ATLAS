
MODULE Fit
	!USE ?
	IMPLICIT NONE

CONTAINS

	SUBROUTINE FitPspec (myid, lon, lat, spec, ntheta, nk, nnu, delta_k, delta_nu, fname_guess, verbose, outdir)
		IMPLICIT NONE

		INTEGER :: ntheta, nk, nnu, myid
		DOUBLE PRECISION :: delta_k, delta_nu
		REAL :: lon, lat
		DOUBLE PRECISION :: spec(ntheta,nk,nnu)
		LOGICAL :: verbose
		CHARACTER(LEN=*) :: fname_guess, outdir
		CHARACTER(LEN=400) :: fname
		INTEGER :: verb, nmodes, ii
		REAL :: dataout(4000)
		INTEGER :: n(1000), k(1000), kstart, kend

		verb = 0
		IF (verbose) verb = 1
		
		WRITE(fname,'(A,A6,I0.5)') &
			TRIM(outdir),"/fits_",myid

!$OMP PARALLEL DO PRIVATE(kstart,kend,ii,nmodes,dataout,n,k)
		DO kstart=4,50,4
			kend = kstart+3
			! fit_wrapper is C code, look for fit_wrapper.c
			! this is the end of FORTRAN code, besides the little bit below that
			! outputs the fitting result to the disk
			CALL fit_wrapper (spec, ntheta, nk, nnu, delta_k, delta_nu, TRIM(fname_guess), verb, nmodes, dataout, n, k, kstart, kend)
!$OMP CRITICAL
      if (verbose) PRINT*, "Modes to write to disk: ",nmodes
			OPEN(4,FILE=TRIM(fname),ACCESS="APPEND",ACTION="WRITE")
			DO ii=1,nmodes
				WRITE(4,'(2F12.6,2I5,4E15.7)') lon, lat, n(ii), k(ii), dataout(ii*4:ii*4+3)
			ENDDO
			CLOSE(4)
!$OMP END CRITICAL
		ENDDO
!$OMP END PARALLEL DO

	END SUBROUTINE FitPspec


END MODULE Fit
